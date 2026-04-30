"""Aggregate per-PDB context tsvs and compute the signed angle
between adjacent C-alpha -> side-chain centroid vectors for each
HHH tripeptide centred on the target amino acid.

For every tripeptide the angle is taken between the LEFT and CENTRE
residue: v_left = CA_L -> centroid_L, v_centre = CA_C -> centroid_C,
about the axis CA_L -> CA_C.

Writes:
    final/angles.tsv      one row per tripeptide
    final/valid_pdbs.txt  PDB IDs that contributed at least one row
"""

import gzip
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np


SIZE_CLASS = {
    "G": "Tiny",  "A": "Tiny",
    "V": "Small", "P": "Small", "S": "Small", "T": "Small", "C": "Small",
    "D": "Intermediate", "L": "Intermediate", "I": "Intermediate", "N": "Intermediate",
    "K": "Large", "E": "Large", "M": "Large", "Q": "Large", "H": "Large",
    "R": "Bulky", "F": "Bulky", "Y": "Bulky", "W": "Bulky",
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

BACKBONE = {"N", "CA", "C", "O", "OXT"}

OUT_HEADER = [
    "pdb_id", "chain", "resseq",
    "left_aa", "center_aa", "right_aa",
    "left_size_class", "angle",
]


def signed_angle_3d(v1, v2, axis, degrees=True):
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    axis = np.array(axis, dtype=float)
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    axis /= np.linalg.norm(axis)
    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)
    angle = np.arctan2(np.dot(cross, axis), dot)
    return np.degrees(angle) if degrees else angle


def parse_pdb_coords(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    out = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            atom = line[12:16].strip()
            chain = line[21]
            resseq = line[22:26].strip()
            try:
                xyz = np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ])
            except ValueError:
                continue
            entry = out.setdefault((chain, resseq), {"ca": None, "sc": []})
            if atom == "CA":
                entry["ca"] = xyz
            elif atom not in BACKBONE:
                entry["sc"].append(xyz)
    return out


def centroid(entry):
    if entry["sc"]:
        return np.mean(entry["sc"], axis=0)
    return entry["ca"]


def read_context(path):
    """Return dict[tag] -> [(role, resname, chain, resseq), ...]"""
    groups = defaultdict(list)
    pdb_id = None
    try:
        fh = open(path)
    except OSError:
        return pdb_id, groups
    with fh:
        header = fh.readline()
        if not header:
            return pdb_id, groups
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            pdb_id = parts[0]
            tag    = parts[1]
            role   = parts[2]
            resn   = parts[3]
            chain  = parts[4]
            resseq = parts[5]
            groups[tag].append((role, resn, chain, resseq))
    return pdb_id, groups


def angles_for_one(args):
    """Worker: one context tsv -> list of output rows."""
    ctx_path, pdb_dir = args
    pdb_id, groups = read_context(ctx_path)
    if pdb_id is None or not groups:
        return []

    pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb.gz")
    if not os.path.exists(pdb_file):
        pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")
    try:
        coords = parse_pdb_coords(pdb_file)
    except (FileNotFoundError, OSError):
        return []

    rows = []
    for tag, members in groups.items():
        by_role = {role: (resn, ch, rs) for (role, resn, ch, rs) in members}
        if not all(k in by_role for k in ("L", "C", "R")):
            continue
        L = by_role["L"]; C = by_role["C"]; R = by_role["R"]

        try:
            el = coords[(L[1], L[2])]
            ec = coords[(C[1], C[2])]
        except KeyError:
            continue
        if el["ca"] is None or ec["ca"] is None:
            continue

        v_left   = centroid(el) - el["ca"]
        v_centre = centroid(ec) - ec["ca"]
        if np.linalg.norm(v_left) < 1e-6 or np.linalg.norm(v_centre) < 1e-6:
            continue
        axis = ec["ca"] - el["ca"]
        if np.linalg.norm(axis) < 1e-6:
            continue

        ang = signed_angle_3d(v_left, v_centre, axis)

        l1 = THREE_TO_ONE.get(L[0]); c1 = THREE_TO_ONE.get(C[0]); r1 = THREE_TO_ONE.get(R[0])
        if l1 is None or r1 is None:
            continue
        size_cls = SIZE_CLASS.get(l1)
        if size_cls is None:
            continue

        rows.append([
            pdb_id, C[1], C[2],
            l1, c1, r1, size_cls, f"{ang:.4f}",
        ])
    return rows


def main():
    if len(sys.argv) < 6:
        sys.exit(
            "usage: compute_angles.py <contexts_dir> <out_angles_tsv> "
            "<target_aa> <pdb_dir> <out_valid_pdbs>"
        )
    ctx_dir, out_angles, target_aa, pdb_dir, out_valid = sys.argv[1:6]
    workers = int(os.environ.get("BET_WORKERS", "8"))

    target_aa = target_aa.upper()
    files = sorted(
        os.path.join(ctx_dir, f) for f in os.listdir(ctx_dir)
        if f.startswith(f"context_for_{target_aa}_in_") and f.endswith(".tsv")
    )

    os.makedirs(os.path.dirname(out_angles) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(out_valid) or ".", exist_ok=True)

    valid = []
    total = 0

    with open(out_angles, "w") as fo:
        fo.write("\t".join(OUT_HEADER) + "\n")
        with Pool(workers) as pool:
            for rows in pool.imap_unordered(
                angles_for_one,
                ((f, pdb_dir) for f in files),
                chunksize=64,
            ):
                if not rows:
                    continue
                pid = rows[0][0]
                valid.append(pid)
                total += len(rows)
                for r in rows:
                    fo.write("\t".join(map(str, r)) + "\n")

    valid.sort()
    with open(out_valid, "w") as fo:
        for pid in valid:
            fo.write(pid + "\n")

    print(f"angle rows: {total}, valid pdbs: {len(valid)}")


if __name__ == "__main__":
    main()
