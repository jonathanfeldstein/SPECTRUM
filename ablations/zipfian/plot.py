#!/usr/bin/env python3
"""Plots ablations/zipfian/rules.csv: precision (symmetry- and prior-corrected) against recall for every
rule SPECTRUM can learn (reference run), one log-log panel per dataset. Red crosses were missed in at least
one default-hyperparameter run. The ring marks the Family rule discussed in the paper. The highest-precision
rule in the bottom recall quartile of each dataset is printed.

Sized for a full-width figure: at 100% zoom on a 7-inch text width, labels print at about 9-10 pt.

usage: python3 ablations/zipfian/plot.py   (needs matplotlib)
Writes ablations/zipfian/figure.png and figure.pdf.
"""
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
NAMES = {"family": "Family", "umls": "UMLS", "kinship": "Kinship", "wn18rr": "WN18RR", "fb15k237": "FB15K-237"}
RINGED = {"family": "mother(A,B) & father(C,B) >> wife(A,C)"}  # the example discussed in the paper
MINED, MISSED = "#6a9bd1", "#e05a5b"

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11, "legend.fontsize": 10,
                     "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
                     "axes.spines.top": False, "axes.spines.right": False})


def y_ticks(ax):
    """Plain-number labels at 1, 2, 5, 10, ... (or decades only on axes spanning more than 2.5 decades)."""
    lo, hi = ax.get_ylim()
    subs = (1,) if math.log10(hi / lo) > 2.5 else (1, 2, 5)
    ticks = [s * 10**e for e in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1)
             for s in subs if lo <= s * 10**e <= hi]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(NullLocator())


def main():
    with open(HERE / "rules.csv") as f:
        rows = list(csv.DictReader(f))
    datasets = [d for d in NAMES if any(r["dataset"] == d for r in rows)]
    if not datasets:
        sys.exit("rules.csv has no rows; run ablations/zipfian/run.py first")

    fig, axes = plt.subplots(1, len(datasets), figsize=(2.3 * len(datasets), 2.6), squeeze=False)
    for ax, d in zip(axes[0], datasets):
        rs = [r for r in rows if r["dataset"] == d]
        ax.set_xscale("log")
        ax.set_yscale("log")
        for missed in (False, True):
            sel = [r for r in rs if (int(r["mined_runs"]) < int(r["repeats"])) == missed]
            ax.scatter([float(r["recall"]) for r in sel], [float(r["fidelity"]) for r in sel],
                       s=18 if missed else 8, marker="x" if missed else "o", color=MISSED if missed else MINED,
                       alpha=1 if missed else 0.6, linewidths=1.2 if missed else 0, zorder=3 if missed else 2)
        n_missed = sum(int(r["mined_runs"]) < int(r["repeats"]) for r in rs)
        # counts in the top-left corner, above the data: raise the top of the axis by a third of its log range
        y_ticks(ax)  # before the headroom, so only the data range is labelled
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi * (hi / lo) ** 0.33)
        ax.text(0.04, 0.97, f"{len(rs) - n_missed:,} mined", transform=ax.transAxes, va="top",
                color="#3a6ea5", fontsize=9.5)
        ax.text(0.04, 0.85, f"{n_missed:,} missed", transform=ax.transAxes, va="top", color=MISSED, fontsize=9.5)
        for r in rs:
            if r["rule"] == RINGED.get(d):
                ax.scatter([float(r["recall"])], [float(r["fidelity"])], s=120, facecolors="none",
                           edgecolors="black", linewidths=1.2, zorder=4)
        tail = [r for r in rs if float(r["recall_pct"]) <= 25]
        if tail:
            best = max(tail, key=lambda r: float(r["fidelity"]))
            print(f"{NAMES[d]:9s} highest precision in bottom recall quartile: {best['rule']}  "
                  f"(recall pct {best['recall_pct']}, precision pct {best['fidelity_pct']}, "
                  f"mined {best['mined_runs']}/{best['repeats']}, returned {best['returned_runs']}/{best['repeats']})")
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_title(NAMES[d])
    axes[0][0].set_ylabel("Precision")
    fig.supxlabel("Recall", fontsize=11, y=0.1)
    handles = [Line2D([], [], ls="", marker="o", ms=5, mfc=MINED, mec="none", label="mined in every run"),
               Line2D([], [], ls="", marker="x", ms=6, mew=1.5, color=MISSED, label="missed in at least one run")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.08),
               handletextpad=0.2, columnspacing=1.5)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"figure.{ext}", dpi=300, bbox_inches="tight")
    print(f"Wrote {HERE / 'figure.png'} and {HERE / 'figure.pdf'}")


if __name__ == "__main__":
    main()
