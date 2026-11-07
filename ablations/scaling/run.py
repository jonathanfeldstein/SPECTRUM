#!/usr/bin/env python3
"""Scaling: SPECTRUM's mining time vs dataset size.

Mines every dataset --repeats times with its own hyperparameters (configs/hyperparameters.tsv,
M = 20 x number of predicates) and records the wall time, M and the number of facts (atoms in
train.db) in results.csv. Runs already there are skipped, so the script can be stopped and resumed.
Run it on an otherwise idle machine: the runtimes are the result.

usage: python3 ablations/scaling/run.py [--datasets ...] [--repeats 5]
Then: python3 ablations/scaling/plot.py
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import ROOT, ResultsWriter, check_binaries, heuristic_M, hyperparameters, read_results, run_miner  # noqa: E402

RESULTS = HERE / "results.csv"
FIELDS = ["dataset", "repeat", "facts", "M", "epsilon", "fragments", "rules", "runtime_s", "status"]
DATASETS = ["umls", "kinship", "family", "wn18rr", "cad", "fb15k237", "yelp"]


def facts(dataset):
    """Number of atoms in train.db (comment lines excluded)."""
    with open(ROOT / "datasets" / dataset / "train.db") as f:
        return sum(1 for line in f if "(" in line and not line.startswith("//"))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", nargs="+", default=DATASETS)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--timeout", type=float, default=7200, help="seconds per run")
    args = parser.parse_args()

    check_binaries(["SPECTRUM"])
    hp = hyperparameters()
    done = {(r["dataset"], int(r["repeat"])) for r in read_results(RESULTS)}
    out = ResultsWriter(RESULTS, FIELDS)
    for dataset in args.datasets:
        n_facts, M, eps = facts(dataset), heuristic_M(dataset), hp[dataset]["epsilon"]
        for repeat in range(1, args.repeats + 1):
            if (dataset, repeat) in done:
                continue
            t0 = time.time()
            res = run_miner("SPECTRUM", dataset, M, eps, hp[dataset], args.timeout)
            res.pop("random_walks", None)
            row = dict(dataset=dataset, repeat=repeat, facts=n_facts, M=M, epsilon=eps, **res)
            out.write(row)
            print(f"{dataset:9s} rep={repeat} facts={n_facts} M={M} runtime={row.get('runtime_s', '')}s "
                  f"[{row['status']}, {time.time() - t0:.0f}s]", flush=True)
    out.close()


if __name__ == "__main__":
    main()
