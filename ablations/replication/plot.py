#!/usr/bin/env python3
"""Plots ablations/replication/results.csv as two log-log panels (mean +/- std over repeats):

  left:   fixed N, runtime vs number of facts (k copies); each dataset's fitted slope is listed
  right:  auto N, the same with epsilon fixed (N shrinks as 1/k)

usage: python3 ablations/replication/plot.py   (needs matplotlib and numpy)
Writes ablations/replication/figure.png and figure.pdf and prints the fitted slopes. Runtime vs
pattern locations (all runs: slope 0.98) is no longer plotted; locations remain in results.csv.
"""
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import read_results  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11, "legend.fontsize": 10,
                     "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "lines.linewidth": 1.6,
                     "lines.markersize": 4, "axes.spines.top": False, "axes.spines.right": False})

DATASET_STYLES = {
    "umls": dict(color="#1f77b4", linestyle="-", label="UMLS"),
    "kinship": dict(color="#ff7f0e", linestyle="--", label="Kinship"),
    "family": dict(color="#2ca02c", linestyle="-.", label="Family"),
    "wn18rr": dict(color="#8c564b", linestyle=(0, (5, 2)), label="WN18RR"),
    "fb15k237": dict(color="#9467bd", linestyle=":", label="FB15K-237"),
}
SERIES = [("fixed_N", r"Fixed $N$  ($\varepsilon_k = \varepsilon/\sqrt{k}$)"),
          ("auto_N", r"Fixed $\varepsilon$  ($N \propto 1/k$)")]


def fit_slope(x, y):
    """Least-squares slope of log10(y) on log10(x) and its standard error."""
    (a, _), cov = np.polyfit(np.log10(x), np.log10(y), 1, cov=True)
    return a, np.sqrt(cov[0, 0])


def plain_ticks(axis):
    """Plain-number labels at powers of ten (0.1, 1, 10)."""
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_locator(NullLocator())


def slope_block(ax, entries):
    """Fitted slopes as a small coloured list in the panel's lower-right corner, which the data leave empty."""
    ax.text(0.97, 0.05 + 0.085 * len(entries), "slope", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9.5, color="0.3")
    for i, (text, color) in enumerate(reversed(entries)):
        ax.text(0.97, 0.05 + 0.085 * i, text, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=9.5, color=color)


def main():
    runs = defaultdict(list)  # (dataset, series, copies) -> rows
    for r in read_results(HERE / "results.csv"):
        if r["status"] == "ok" and r["runtime_s"]:
            runs[(r["dataset"], r["series"], int(r["copies"]))].append(r)
    datasets = [d for d in DATASET_STYLES if any(k[0] == d for k in runs)]

    fig, axes = plt.subplots(1, 2, figsize=(7.7, 3.3))
    axes[1].sharey(axes[0])
    for ax, (series, title) in zip(axes, SERIES):
        slopes = []
        for d in datasets:
            ks = sorted(k for (dd, s, k) in runs if dd == d and s == series)
            if not ks:
                continue
            x = np.array([float(runs[(d, series, k)][0]["facts"]) for k in ks])
            t = [np.array([float(r["runtime_s"]) for r in runs[(d, series, k)]]) for k in ks]
            y, yerr = np.array([v.mean() for v in t]), np.array([v.std() for v in t])
            style = DATASET_STYLES[d]
            if len(ks) >= 3:
                a, a_err = fit_slope(x, y)
                slopes.append((f"{style['label']}  {a:.2f}", style["color"]))
                print(f"{series:7s} {d:8s} slope = {a:.3f} ± {a_err:.3f}")
            ax.errorbar(x, y, yerr=yerr, marker="o", capsize=0, elinewidth=1, **style)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_minor_locator(NullLocator())
        plain_ticks(ax.yaxis)
        ax.set_xlabel("Facts")
        ax.set_title(title)
        slope_block(ax, slopes)
    axes[0].set_ylabel("Runtime (s)")

    handles = [Line2D([], [], marker="o", **DATASET_STYLES[d]) for d in datasets]
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 1.08),
               handlelength=2.4, columnspacing=1.5)
    fig.tight_layout(w_pad=2)
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"figure.{ext}", dpi=300, bbox_inches="tight")
    print(f"Wrote {HERE / 'figure.png'} and {HERE / 'figure.pdf'}")


if __name__ == "__main__":
    main()
