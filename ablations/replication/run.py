#!/usr/bin/env python3
"""Replication: SPECTRUM's runtime as one dataset is copied k times (k disjoint copies of its graph).

Copying multiplies the number of nodes, facts and every local neighbourhood count by exactly k while
keeping each copy's structure, so it isolates the effect of dataset size. Two series per dataset:

  fixed_N: epsilon_k = epsilon / sqrt(k). Since N = ceil(M ln M / (epsilon^2 |V|)), this keeps the
           number of paths per node N the same for every k; the mining analysis then predicts
           runtime proportional to k.
  auto_N:  epsilon fixed, so N shrinks as 1/k and the total path budget |V| N ~ M ln M / epsilon^2
           stays constant; runtime is predicted to stay flat while that budget binds (dense graphs)
           and to grow with k where nodes' neighbourhoods are explored exhaustively anyway.

M, D, M_multi and min_fidelity are the dataset's values (configs/hyperparameters.tsv; M = 20 x number of
predicates, which copying does not change). Copy 1 keeps the original constant names, copy i > 1 appends
"__c<i>", so k = 1 is the original train.db. Copies are written to ablations/replication/data/.
Copying is not exactly neutral for the rules: a pattern becomes a rule template only with more than 5
locations (CandidateRules.cpp, get_constrained_fragments), so patterns with 3-5 matches per copy pass
from k = 2 on and the candidate rule count grows slightly until k >= 6. Mining is unaffected.
Runs single-threaded by default so the runtimes are not dominated by thread start-up on small graphs.
Runs already in results.csv are skipped, so the script can be stopped and resumed.

usage: python3 ablations/replication/run.py [--datasets ...] [--copies 1 2 4 8 16] [--repeats 3] [--threads 1]
Then: python3 ablations/replication/plot.py
"""
import argparse
import math
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import ROOT, ResultsWriter, check_binaries, heuristic_M, hyperparameters, read_results, run_miner  # noqa: E402

RESULTS = HERE / "results.csv"
FIELDS = ["dataset", "series", "copies", "repeat", "facts", "M", "epsilon", "random_walks", "locations",
          "fragments", "rules", "mining_s", "runtime_s", "status"]
DATASETS = ["family", "umls", "kinship", "wn18rr"]


def replicate(dataset, k):
    """Writes k disjoint copies of datasets/<dataset>/train.db (once) and returns (folder, facts)."""
    folder = HERE / "data" / f"{dataset}_x{k}"
    facts = [l.strip() for l in (ROOT / "datasets" / dataset / "train.db").read_text().splitlines()
             if "(" in l and not l.startswith("//")]
    target = folder / "train.db"
    if not target.exists():
        folder.mkdir(parents=True, exist_ok=True)
        tmp = target.with_suffix(".tmp")
        with open(tmp, "w") as f:
            for i in range(1, k + 1):
                suffix = "" if i == 1 else f"__c{i}"
                for fact in facts:
                    pred, args = fact[:fact.index("(")], fact[fact.index("(") + 1:fact.rindex(")")]
                    f.write(f"{pred}({', '.join(a.strip() + suffix for a in args.split(','))})\n")
        tmp.rename(target)
    return folder, k * len(facts)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", nargs="+", default=DATASETS)
    parser.add_argument("--copies", nargs="+", type=int, default=[1, 2, 4, 8, 16])
    parser.add_argument("--series", nargs="+", default=["fixed_N", "auto_N"], choices=["fixed_N", "auto_N"])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--threads", type=int, default=1, help="OMP_NUM_THREADS for every run")
    parser.add_argument("--timeout", type=float, default=7200, help="seconds per run")
    args = parser.parse_args()

    check_binaries(["SPECTRUM"])
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    hp = hyperparameters()
    done = {(r["dataset"], r["series"], int(r["copies"]), int(r["repeat"])) for r in read_results(RESULTS)}
    out = ResultsWriter(RESULTS, FIELDS)
    for dataset in args.datasets:
        M = heuristic_M(dataset)
        for k in args.copies:
            folder, n_facts = replicate(dataset, k)
            for series in args.series:
                eps = hp[dataset]["epsilon"] / math.sqrt(k) if series == "fixed_N" else hp[dataset]["epsilon"]
                for repeat in range(1, args.repeats + 1):
                    if (dataset, series, k, repeat) in done:
                        continue
                    t0 = time.time()
                    res = run_miner("SPECTRUM", dataset, M, f"{eps:.10g}", hp[dataset], args.timeout,
                                    data_folder=folder, details=True)
                    row = dict(dataset=dataset, series=series, copies=k, repeat=repeat, facts=n_facts, M=M,
                               epsilon=f"{eps:.6g}", **res)
                    out.write(row)
                    print(f"{dataset:8s} {series:7s} k={k:<3d} rep={repeat} N={row.get('random_walks', '')} "
                          f"locations={row.get('locations', '')} runtime={row.get('runtime_s', '')}s "
                          f"[{row['status']}, {time.time() - t0:.0f}s]", flush=True)
    out.close()


if __name__ == "__main__":
    main()
