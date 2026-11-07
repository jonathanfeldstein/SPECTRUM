#!/usr/bin/env python3
"""Plots ablations/scaling/results.csv: SPECTRUM's mining runtime against dataset size
(log-log, mean +/- std over repeats) with a least-squares power-law fit y = c * x^a.

usage: python3 ablations/scaling/plot.py   (needs matplotlib and numpy)
Writes ablations/scaling/figure.png and figure.pdf.
"""
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 12.5})  # 25% above the default, for the figure's size in the paper
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import read_results  # noqa: E402

NAMES = {"umls": "UMLS", "kinship": "Kinship", "family": "Family", "wn18rr": "WN18RR",
         "cad": "CAD", "fb15k237": "FB15K-237", "yelp": "Yelp"}
# label offset in points (dx, dy) and alignment, chosen so labels clear the points and the fit line
LABELS = {"umls": (0, 14, "center"), "family": (9, 0, "left"), "wn18rr": (9, 0, "left"),
          "cad": (-9, 0, "right"), "fb15k237": (-9, 0, "right"), "yelp": (0, -15, "center")}
# labels with no room beside their point: placed in empty space (data coordinates) with a leader line
LEADER_LABELS = {"kinship": (2.6e3, 3.0)}


def load():
    per = defaultdict(lambda: dict(facts=None, M=None, t=[]))
    for r in read_results(HERE / "results.csv"):
        if r["status"] != "ok" or not r["runtime_s"]:
            continue
        d = per[r["dataset"]]
        d["facts"], d["M"] = int(r["facts"]), int(r["M"])
        d["t"].append(float(r["runtime_s"]))
    return per


def main():
    per = load()
    datasets = sorted(per, key=lambda d: per[d]["facts"])
    x = np.array([per[d]["facts"] for d in datasets], dtype=float)
    t = [np.array(per[d]["t"]) for d in datasets]
    y, yerr = np.array([v.mean() for v in t]), np.array([v.std() for v in t])

    # fit log10(y) = a log10(x) + b; slope uncertainty = standard error from the fit covariance
    (a, b), cov = np.polyfit(np.log10(x), np.log10(y), 1, cov=True)
    a_err = np.sqrt(cov[0, 0])

    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    xs = np.logspace(np.log10(x.min()) - 0.1, np.log10(x.max()) + 0.1, 50)
    ax.plot(xs, 10 ** b * xs ** a, "k--", linewidth=2.5, dashes=(5, 3),
            label=rf"fit: $y \propto x^{{{a:.2f} \pm {a_err:.2f}}}$", zorder=1)
    ax.errorbar(x, y, yerr=yerr, fmt="o", color="#d6403a", markersize=10, capsize=0, zorder=2)
    for d, xi, yi in zip(datasets, x, y):
        dx, dy, ha = LABELS.get(d, (7, 5, "left"))
        if d in LEADER_LABELS:
            ax.annotate(NAMES[d], (xi, yi), xytext=LEADER_LABELS[d], textcoords="data", ha="left", va="center",
                        fontsize=12.5, arrowprops=dict(arrowstyle="-", color="0.35", linewidth=0.8, shrinkB=6))
            continue
        ax.annotate(NAMES.get(d, d), (xi, yi), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va="center", fontsize=12.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Dataset size (number of facts)")
    ax.set_ylabel("Mining runtime (s)")
    ax.set_xlim(xs[0] / 1.3, xs[-1] * 1.3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(loc="lower right", bbox_to_anchor=(1.03, 0.08), frameon=False, handlelength=2.6)  # two whole dashes; empty area below WN18RR
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"figure.{ext}", dpi=300)
    print(f"slope a = {a:.3f} ± {a_err:.3f}")
    print(f"Wrote {HERE / 'figure.png'} and {HERE / 'figure.pdf'}")


if __name__ == "__main__":
    main()
