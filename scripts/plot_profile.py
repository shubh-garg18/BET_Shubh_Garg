"""Density plot of side-chain orientation angles in HHH XRX
tripeptides centred on Arginine, grouped by the size class of the
left neighbour.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


CLASS_ORDER = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]
CLASS_COLOR = {
    "Tiny":         "#dcd5c0",
    "Small":        "#f5c46d",
    "Intermediate": "#ef8a4a",
    "Large":        "#d94c2a",
    "Bulky":        "#9c1024",
}


def wrap_to_180(x):
    return ((x + 180.0) % 360.0) - 180.0


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: make_plot.py <angles.tsv> <out.png> [valid_pdbs.txt]")
    angles_tsv = sys.argv[1]
    out_png    = sys.argv[2]
    valid_path = sys.argv[3] if len(sys.argv) > 3 else None

    df = pd.read_csv(angles_tsv, sep="\t")
    df = df[df["left_size_class"].isin(CLASS_ORDER)].copy()
    df["angle"] = wrap_to_180(df["angle"].astype(float))

    n_angles = len(df)
    n_pdbs = 0
    if valid_path and os.path.exists(valid_path):
        with open(valid_path) as fh:
            n_pdbs = sum(1 for ln in fh if ln.strip())

    grid = np.linspace(-180, 180, 1000)

    fig, ax = plt.subplots(figsize=(10, 6))
    bg = "#bdbdbd"
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)

    for cls in CLASS_ORDER:
        vals = df.loc[df["left_size_class"] == cls, "angle"].to_numpy()
        if len(vals) < 2:
            continue
        kde = gaussian_kde(vals, bw_method="silverman")
        ax.plot(grid, kde(grid), color=CLASS_COLOR[cls], linewidth=1.6, label=cls)

    for x in np.arange(-180, 181, 50):
        ax.axvline(x, color="blue", linestyle=":", linewidth=0.5, alpha=0.7)

    ax.set_xlim(-180, 180)
    ax.set_xticks(np.arange(-180, 181, 50))
    ax.set_xlabel(r"Angle between adjacent C-$\alpha$ $\to$ Centroid vectors [°]")
    ax.set_ylabel("Norm. Freq. [A.U.]")
    title = f"Tripeptide (XRX) in Helix  (n = {n_angles}"
    if n_pdbs:
        title += f", pdbs = {n_pdbs}"
    title += ")"
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    leg = ax.legend(title="Left-neighbour size", loc="upper left",
                    facecolor="white", edgecolor="black")
    leg.get_frame().set_linewidth(0.8)

    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, facecolor=fig.get_facecolor())
    print(f"saved {out_png}  (n={n_angles}, pdbs={n_pdbs})")


if __name__ == "__main__":
    main()
