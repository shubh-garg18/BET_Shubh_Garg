# BET-105

Snakemake pipeline for the assignment. Takes a folder of `.pdb.gz`
files, runs STRIDE on each, picks out every **HHH XRX** tripeptide
where the centre residue is arginine, and computes the signed angle
between adjacent C-α → side-chain centroid vectors (the left and the
centre residue of the tripeptide, about the CA→CA axis). Plots the
result as a KDE grouped by the size class of the left neighbour.

## Plot

![ss profile](final/ss_profile_HHH_for_arg_with_valid_runs.png)

PDB IDs that contributed to the plot:
[`final/valid_pdbs.txt`](final/valid_pdbs.txt)

## Run

Place `.pdb.gz` files inside a `pdbs/` folder at the repo root, then:

```bash
snakemake --cores 8
```

The plot is written to
`final/ss_profile_HHH_for_arg_with_valid_runs.png`.

The conda environment with snakemake, stride and the python deps is in
`environment.yml`:

```bash
conda env create -f environment.yml
conda activate bet-Shubh_Garg
```

## Layout

```
Snakefile                5 rules: unzip_pdb -> run_stride
                         -> extract_context -> calculate_angles -> plot_angles
config.yaml              default target_aa, ss_pattern, paths
scripts/
  extract_contexts.py    per-PDB: stride ASG -> tripeptide context tsv
  calc_angles.py         aggregate contexts + read PDB coords -> angles.tsv,
                         valid_pdbs.txt (multiprocessing)
  plot_profile.py        KDE plot grouped by left-neighbour size
final/
  angles.tsv
  valid_pdbs.txt
  ss_profile_HHH_for_arg_with_valid_runs.png
```

## Size classes (left neighbour)

| Class        | Residues          |
|--------------|-------------------|
| Tiny         | G, A              |
| Small        | V, P, S, T, C     |
| Intermediate | D, L, I, N        |
| Large        | K, E, M, Q, H     |
| Bulky        | R, F, Y, W        |
