configfile: "config.yaml"

import glob
import os

PDB_DIR    = config.get("pdb_dir", "pdbs")
TARGET_AA  = config.get("target_aa", "ARG").upper()
SS_PATTERN = config.get("ss_pattern", "HHH").upper()
STRIDE_BIN = config.get("stride_bin", "stride")
LIST_FILE  = config.get("list_file", "")


def _strip_ext(name):
    for suf in (".pdb.gz", ".pdb"):
        if name.endswith(suf):
            return name[: -len(suf)]
    return name


if LIST_FILE and os.path.exists(LIST_FILE):
    with open(LIST_FILE) as fh:
        PDB_IDS = [_strip_ext(ln.strip()) for ln in fh if ln.strip()]
else:
    PDB_IDS = sorted(
        _strip_ext(os.path.basename(p))
        for p in glob.glob(os.path.join(PDB_DIR, "*.pdb.gz"))
    )


PLOT_PNG   = "results/ss_profile_HHH_for_arg_with_valid_runs.png"
ANGLES_TSV = "results/angles.tsv"
VALID_PDBS = "results/valid_pdbs.txt"


rule all:
    input:
        PLOT_PNG,
        ANGLES_TSV,
        VALID_PDBS,


rule unzip_pdb:
    input:
        pdb_gz = lambda w: f"{PDB_DIR}/{w.pdb}.pdb.gz",
    output:
        pdb = temp("unzipped_pdbs/{pdb}.pdb"),
    shell:
        "zcat {input.pdb_gz} > {output.pdb}"


rule run_stride:
    input:
        pdb = "unzipped_pdbs/{pdb}.pdb",
    output:
        ss = "stride-out/{pdb}.ss.txt",
    params:
        stride = STRIDE_BIN,
    shell:
        '{params.stride} {input.pdb} > {output.ss} 2>/dev/null '
        '|| echo "FAILED {wildcards.pdb}" > {output.ss}'


rule extract_context:
    input:
        ss = "stride-out/{pdb}.ss.txt",
    output:
        ctx = "contexts/context_for_{aa}_in_{pdb}.tsv",
    shell:
        "python scripts/extract_contexts.py {input.ss} {output.ctx} {wildcards.aa}"


rule calculate_angles:
    input:
        contexts = expand(
            "contexts/context_for_{aa}_in_{pdb}.tsv",
            pdb=PDB_IDS, aa=[TARGET_AA],
        ),
    output:
        angles     = ANGLES_TSV,
        valid_pdbs = VALID_PDBS,
    params:
        aa      = TARGET_AA,
        pdb_dir = PDB_DIR,
    shell:
        "python scripts/calc_angles.py contexts {output.angles} "
        "{params.aa} {params.pdb_dir} {output.valid_pdbs}"


rule plot_angles:
    input:
        angles     = ANGLES_TSV,
        valid_pdbs = VALID_PDBS,
    output:
        png = PLOT_PNG,
    shell:
        "python scripts/plot_profile.py {input.angles} {output.png} {input.valid_pdbs}"
