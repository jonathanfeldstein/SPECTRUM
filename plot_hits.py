#!/usr/bin/env python3
"""Plots Hit@1..Hit@10 for every dataset evaluated by ./reproduce.sh.

Reads Hit@1..10 from results/<dataset>/eval.log and writes
results/hits_at_k.png and results/hits_at_k.pdf.

usage: python3 plot_hits.py   (run ./reproduce.sh first; needs matplotlib)
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = Path(__file__).resolve().parent / "results"
# Sized for a half-width figure: at 100% zoom on a 3.4-inch column, labels print at about 9 pt.
plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "legend.fontsize": 9.5,
                     "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
                     "axes.spines.top": False, "axes.spines.right": False})
STYLES = {
    "family": dict(label="Family", color="#4f9a3a", marker="^", linestyle="-."),
    "umls": dict(label="UMLS", color="#3b75af", marker="o", linestyle="-"),
    "kinship": dict(label="Kinship", color="#e8833a", marker="s", linestyle="--"),
    "wn18rr": dict(label="WN18RR", color="#8c564b", marker="D", linestyle="--"),
    "fb15k237": dict(label="FB15K-237", color="#9467bd", marker="v", linestyle="-"),
}


def hits(eval_log):
    """Hit@1..10 from an eval.log."""
    return [float(v) for v in re.findall(r"^Hit@\d+:\s+(\S+)", eval_log.read_text(), re.M)][:10]


def main():
    fig, ax = plt.subplots(figsize=(4.6, 3.1))
    ks = list(range(1, 11))
    for dataset, style in STYLES.items():
        log = RESULTS / dataset / "eval.log"
        if not log.exists():
            continue
        ax.plot(ks, hits(log), linewidth=1.6, markersize=5, **style)
    ax.set_xticks(ks)
    ax.set_xlabel("$k$")
    ax.set_ylabel("Hit@$k$")
    ax.set_ylim(0.2, 1.02)
    ax.grid(True, axis="y", color="#e5e5e5", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", frameon=False, ncol=2, columnspacing=1.0, handlelength=2.2)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(RESULTS / f"hits_at_k.{ext}", dpi=300)
    print(f"Wrote {RESULTS / 'hits_at_k.png'} and {RESULTS / 'hits_at_k.pdf'}")


if __name__ == "__main__":
    main()
