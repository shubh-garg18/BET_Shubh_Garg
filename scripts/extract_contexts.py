"""Per-PDB context extractor.

Reads a STRIDE assignment file, walks ASG records, and for every
HHH tripeptide whose centre residue equals the target amino acid,
emits three rows -- one per residue in the window -- so the angle
calculation step can pick coordinates without re-parsing STRIDE.
"""

import os
import sys


HEADER = [
    "pdb_id", "tripeptide_pos",
    "role", "resname", "chain", "resseq", "ord",
    "ss_code", "ss_name", "phi", "psi", "area",
    "tripeptide_seq", "ss_pattern",
]

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_asg(line):
    p = line.split()
    if len(p) < 10 or p[0] != "ASG":
        return None
    return {
        "resname": p[1],
        "chain":   p[2],
        "resseq":  p[3],
        "ord":     p[4],
        "ss_code": p[5],
        "ss_name": p[6],
        "phi":     p[7],
        "psi":     p[8],
        "area":    p[9],
    }


def read_stride(path):
    rows = []
    try:
        fh = open(path)
    except OSError:
        return rows
    with fh:
        for line in fh:
            if line.startswith("ASG"):
                r = parse_asg(line)
                if r is not None:
                    rows.append(r)
    return rows


def main():
    if len(sys.argv) < 4:
        sys.exit("usage: make_contexts.py <stride.ss.txt> <out.tsv> <three_letter_aa>")
    in_ss, out_tsv, target = sys.argv[1], sys.argv[2], sys.argv[3].upper()

    pdb_id = os.path.basename(in_ss).split(".")[0]
    residues = read_stride(in_ss)

    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w") as fo:
        fo.write("\t".join(HEADER) + "\n")

        for i in range(1, len(residues) - 1):
            c = residues[i]
            if c["resname"] != target:
                continue
            l, r = residues[i - 1], residues[i + 1]
            if l["chain"] != c["chain"] or r["chain"] != c["chain"]:
                continue
            ss_pat = l["ss_code"] + c["ss_code"] + r["ss_code"]
            if ss_pat != "HHH":
                continue

            seq = "".join(THREE_TO_ONE.get(x["resname"], "X") for x in (l, c, r))
            tag = f"{pdb_id}:{c['chain']}:{c['resseq']}"

            for role, res in (("L", l), ("C", c), ("R", r)):
                fo.write("\t".join([
                    pdb_id, tag, role,
                    res["resname"], res["chain"], res["resseq"], res["ord"],
                    res["ss_code"], res["ss_name"],
                    res["phi"], res["psi"], res["area"],
                    seq, ss_pat,
                ]) + "\n")


if __name__ == "__main__":
    main()
