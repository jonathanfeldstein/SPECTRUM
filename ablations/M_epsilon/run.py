#!/usr/bin/env python3
"""Ablation: SPECTRUM's fragments found, rules mined and runtime as M and epsilon vary.

M ablation:       epsilon = 0.01, M = 2, 4, 8, ... up to 16 x the heuristic M, plus the heuristic M
Epsilon ablation: M = heuristic M (20 x number of predicates), epsilon in --epsilons
All other hyperparameters are the dataset's values in configs/hyperparameters.tsv. Results are
appended to results.csv; runs already there are skipped, so the script can be stopped and resumed.

usage: python3 ablations/M_epsilon/run.py [--datasets umls kinship family wn18rr] [--repeats 3]
Then: python3 ablations/M_epsilon/plot.py
"""
import argparse
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import ResultsWriter, check_binaries, heuristic_M, hyperparameters, read_results, run_miner  # noqa: E402

RESULTS = HERE / "results.csv"
FIELDS = ["dataset", "M", "epsilon", "repeat", "heuristic_M", "random_walks", "fragments", "rules", "runtime_s", "status"]
BASE_EPSILON = 0.01
EPSILONS = [0.1, 0.05, 0.01, 0.005, 0.001]


def settings(dataset, epsilons):
    """(M, epsilon) pairs for one dataset, cheapest first."""
    Mh = heuristic_M(dataset)
    Ms = [2 ** k for k in range(1, math.ceil(math.log2(16 * Mh)) + 1)] + [Mh]
    pairs = {(M, BASE_EPSILON) for M in Ms} | {(Mh, e) for e in epsilons}
    return sorted(pairs, key=lambda p: p[0] * math.log(p[0]) / p[1] ** 2)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", nargs="+", default=["umls", "kinship", "family", "wn18rr"])
    parser.add_argument("--epsilons", nargs="+", type=float, default=EPSILONS)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=7200, help="seconds per run")
    args = parser.parse_args()

    check_binaries(["SPECTRUM"])
    hp = hyperparameters()
    done = {(r["dataset"], int(r["M"]), float(r["epsilon"]), int(r["repeat"])) for r in read_results(RESULTS)}
    out = ResultsWriter(RESULTS, FIELDS)
    for dataset in args.datasets:
        for M, eps in settings(dataset, args.epsilons):
            for repeat in range(1, args.repeats + 1):
                if (dataset, M, eps, repeat) in done:
                    continue
                t0 = time.time()
                row = dict(dataset=dataset, M=M, epsilon=eps, repeat=repeat, heuristic_M=heuristic_M(dataset),
                           **run_miner("SPECTRUM", dataset, M, eps, hp[dataset], args.timeout))
                out.write(row)
                print(f"{dataset:9s} M={M:<6} eps={eps:<6} rep={repeat} fragments={row.get('fragments', '')} "
                      f"rules={row.get('rules', '')} runtime={row.get('runtime_s', '')}s [{row['status']}, "
                      f"{time.time() - t0:.0f}s]", flush=True)
    out.close()


if __name__ == "__main__":
    main()
