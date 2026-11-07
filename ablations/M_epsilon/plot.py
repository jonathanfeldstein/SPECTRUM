#!/usr/bin/env python3
"""Plots ablations/M_epsilon/results.csv as one row of four panels with a shared legend above:
rules and runtime as M varies (epsilon = 0.01; stars mark the heuristic M), then rules and runtime as
epsilon varies (M = heuristic M; epsilon axis inverted, so the sampling budget grows to the right in
every panel). Mean +/- std over repeats. Fragment counts stay in results.csv but are not plotted.

Sized for a full-width figure: at 100% zoom on a 7-inch text width, labels print at about 9-10 pt.

usage: python3 ablations/M_epsilon/plot.py   (needs matplotlib and numpy)
Writes ablations/M_epsilon/figure.png and figure.pdf.
"""
from collections import defaultdict
from pathlib import Path
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import read_results  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11, "legend.fontsize": 10,
                     "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "lines.linewidth": 1.6,
                     "lines.markersize": 3.5, "axes.spines.top": False, "axes.spines.right": False})

# One fixed colour + line style per dataset
DATASET_STYLES = {
    "umls": dict(color="#1f77b4", linestyle="-", label="UMLS"),
    "kinship": dict(color="#ff7f0e", linestyle="--", label="Kinship"),
    "family": dict(color="#2ca02c", linestyle="-.", label="Family"),
    "wn18rr": dict(color="#8c564b", linestyle=(0, (5, 2)), label="WN18RR"),
    "fb15k237": dict(color="#9467bd", linestyle=":", label="FB15K-237"),
}
BASE_EPSILON = 0.01
COLUMNS = [("rules", "Rules", "Rules"), ("runtime_s", "Runtime (s)", "Runtime")]


def load():
    """(dataset, M, epsilon) -> {metric: [values over repeats]}, and dataset -> heuristic M."""
    data, heuristic = defaultdict(lambda: defaultdict(list)), {}
    for r in read_results(HERE / "results.csv"):
        if r["status"] != "ok":
            continue
        key = (r["dataset"], int(r["M"]), float(r["epsilon"]))
        heuristic[r["dataset"]] = int(r["heuristic_M"])
        for metric, _, _ in COLUMNS:
            if r[metric] != "":
                data[key][metric].append(float(r[metric]))
    return data, heuristic


def series(data, dataset, metric, vary, heuristic):
    """x, mean, std for one dataset: vary='M' (epsilon fixed) or 'epsilon' (M = heuristic)."""
    if vary == "M":
        keys = sorted(k for k in data if k[0] == dataset and k[2] == BASE_EPSILON)
        x = [k[1] for k in keys]
    else:
        keys = sorted(k for k in data if k[0] == dataset and k[1] == heuristic[dataset])
        x = [k[2] for k in keys]
    values = [data[k][metric] for k in keys]
    # zero counts (e.g. no rules at tiny M) cannot be shown on a log axis
    keep = [i for i, v in enumerate(values) if v and np.mean(v) > 0]
    return (np.array(x)[keep], np.array([np.mean(values[i]) for i in keep]),
            np.array([np.std(values[i]) for i in keep]))


def log_ticks(axis, lo, hi, formatter, subs=(1,)):
    """Labels at subs x powers of ten (powers of ten only by default), as plain numbers or via formatter."""
    axis.set_major_locator(FixedLocator([s * 10.0**e for e in range(math.floor(math.log10(lo)),
                                                                    math.ceil(math.log10(hi)) + 1) for s in subs]))
    axis.set_major_formatter(FuncFormatter(formatter))
    axis.set_minor_locator(NullLocator())


def main():
    data, heuristic = load()
    datasets = [d for d in DATASET_STYLES if d in heuristic]
    fig, axes = plt.subplots(1, 4, figsize=(11.5, 2.9))
    axes[2].sharey(axes[0])  # same scales for the M and epsilon panels of each metric
    axes[3].sharey(axes[1])
    panels = [("M", "rules"), ("M", "runtime_s"), ("epsilon", "rules"), ("epsilon", "runtime_s")]
    for ax, (vary, metric) in zip(axes, panels):
        for d in datasets:
            style = DATASET_STYLES[d]
            x, mean, std = series(data, d, metric, vary, heuristic)
            ax.errorbar(x, mean, yerr=std, marker="o", capsize=0, elinewidth=1, **style)
            if vary == "M":  # heuristic M
                star = data.get((d, heuristic[d], BASE_EPSILON), {}).get(metric)
                if star:
                    ax.scatter([heuristic[d]], [np.mean(star)], marker="*", s=130, color=style["color"],
                               edgecolors="black", linewidths=0.6, zorder=5)
        ax.set_yscale("log")
        if vary == "M":
            ax.set_xscale("log", base=2)
            ax.xaxis.set_major_locator(FixedLocator([2.0**e for e in range(1, 18, 4)]))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"$2^{{{int(round(math.log2(v)))}}}$"))
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xlabel("$M$  ($\\varepsilon = 0.01$)")
        else:
            ax.set_xscale("log")
            ax.invert_xaxis()
            ax.set_xticks([0.1, 0.01, 0.001])
            ax.set_xticklabels(["0.1", "0.01", "0.001"])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xlabel("$\\varepsilon$  (heuristic $M$)")
        ax.set_ylabel("Rules" if metric == "rules" else "Runtime (s)")
    for ax, (_, metric) in zip(axes, panels):  # after all panels are drawn, so shared limits are final
        lo, hi = ax.get_ylim()
        log_ticks(ax.yaxis, lo, hi, lambda v, _: f"{v:g}", subs=(1, 2, 5) if metric == "rules" else (1,))
        ax.set_ylim(lo, hi)
    handles = [Line2D([], [], marker="o", **DATASET_STYLES[d]) for d in datasets]
    handles.append(Line2D([], [], ls="", marker="*", ms=11, mfc="white", mec="black", label="heuristic $M$"))
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 1.08),
               handlelength=2.2, columnspacing=1.2)
    fig.tight_layout(w_pad=1.5)
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"figure.{ext}", dpi=300, bbox_inches="tight")
    print(f"Wrote {HERE / 'figure.png'} and {HERE / 'figure.pdf'}")


if __name__ == "__main__":
    main()
