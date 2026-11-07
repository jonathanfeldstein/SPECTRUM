#!/usr/bin/env python3
"""Plots ablations/random_walk_comparison/results.csv: rules mined and runtime vs epsilon, pattern mining (SPECTRUM)
vs random walks, one column per dataset (mean +/- std over repeats). Fragment counts are not plotted:
the random-walk variant runs with force_cyclic off, so its counts include non-cyclic fragments that
SPECTRUM does not count.

Sized for a full-width figure: at 100% zoom on a 7-inch text width, labels print at about 9-10 pt.

usage: python3 ablations/random_walk_comparison/plot.py   (needs matplotlib and numpy)
Writes ablations/random_walk_comparison/figure.png and figure.pdf.
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
DATASETS = [("family", "Family"), ("kinship", "Kinship"), ("umls", "UMLS"), ("wn18rr", "WN18RR"), ("fb15k237", "FB15K-237")]
VARIANTS = [("SPECTRUM", "SPECTRUM", dict(color="#3b75af", marker="o", linestyle="-")),
            ("RW", "Random walks", dict(color="#e8833a", marker="s", linestyle="--"))]
ROWS = [("rules", "Rules"), ("runtime_s", "Runtime (s)")]

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11, "legend.fontsize": 10,
                     "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
                     "axes.spines.top": False, "axes.spines.right": False})


def load():
    data = defaultdict(list)  # (dataset, variant, metric, eps) -> values
    with open(HERE / "results.csv") as f:
        for r in csv.DictReader(f):
            if r["status"] != "ok":
                continue
            for metric, _ in ROWS:
                if r[metric] != "":
                    data[(r["dataset"], r["variant"], metric, float(r["epsilon"]))].append(float(r[metric]))
    return data


def decade_ticks(ax):
    """Plain-number labels at powers of ten on a log y axis."""
    lo, hi = ax.get_ylim()
    ax.yaxis.set_major_locator(FixedLocator([10.0**e for e in range(math.floor(math.log10(lo)),
                                                                     math.ceil(math.log10(hi)) + 1)]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(NullLocator())


def main():
    data = load()
    datasets = [(d, name) for d, name in DATASETS if any(k[0] == d for k in data)]
    fig, axes = plt.subplots(2, len(datasets), figsize=(2.3 * len(datasets), 4.4), sharex=True, squeeze=False)
    for col, (dataset, name) in enumerate(datasets):
        for row, (metric, label) in enumerate(ROWS):
            ax = axes[row][col]
            for variant, vlabel, style in VARIANTS:
                eps = sorted(e for (d, v, m, e) in data if (d, v, m) == (dataset, variant, metric))
                mean = np.array([np.mean(data[(dataset, variant, metric, e)]) for e in eps])
                std = np.array([np.std(data[(dataset, variant, metric, e)]) for e in eps])
                ax.errorbar(eps, mean, yerr=std, label=vlabel, capsize=2, markersize=4, linewidth=1.4, **style)
            ax.set_xscale("log")
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xticks([0.01, 0.1])
            ax.set_xticklabels(["0.01", "0.1"])
            if metric == "rules":
                ax.set_ylim(bottom=0)  # random walks find no rules at large epsilon on small datasets
                ax.set_title(name)
            else:
                ax.set_yscale("log")
                decade_ticks(ax)
            if col == 0:
                ax.set_ylabel(label)
    fig.supxlabel("ε", fontsize=11, y=0.04)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.05))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"figure.{ext}", dpi=300, bbox_inches="tight")
    print(f"Wrote {HERE / 'figure.png'} and {HERE / 'figure.pdf'}")


if __name__ == "__main__":
    main()
