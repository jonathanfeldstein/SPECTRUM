#!/usr/bin/env python3
"""Ablation: SPECTRUM's pattern mining vs random walks, as a function of epsilon.

For every (dataset, epsilon, repeat) this runs build/SPECTRUM (pattern mining) and
build/SPECTRUM_RW (random walks; force_cyclic off) with the dataset's hyperparameters
from configs/hyperparameters.tsv (only epsilon varies) and M = 20 x number of predicates in
the .info file, and records the number of fragments (patterns) found, candidate rules mined
and the mining wall time in results.csv. Runs already in results.csv are skipped, so the
script can be interrupted and restarted.

usage: python3 ablations/random_walk_comparison/run.py [--datasets family umls ...] [--epsilons 0.5 0.1 ...]
                                                       [--repeats 3] [--variants SPECTRUM RW] [--timeout 7200]
Then: python3 ablations/random_walk_comparison/plot.py
"""
import argparse
import csv
import json
import os
import re
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.csv"
DATASETS = ["family", "kinship", "umls", "wn18rr", "fb15k237"]
EPSILONS = [0.5, 0.1, 0.05, 0.01, 0.005]
BINARIES = {"SPECTRUM": ROOT / "build" / "SPECTRUM", "RW": ROOT / "build" / "SPECTRUM_RW"}
FIELDS = ["dataset", "variant", "epsilon", "repeat", "M", "random_walks", "fragments", "rules", "runtime_s", "status"]


def hyperparameters():
    hp = {}
    for line in (ROOT / "configs" / "hyperparameters.tsv").read_text().splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("dataset"):
            continue
        name, eps, d, m_multi, min_fidelity = line.split("\t")
        hp[name] = dict(D=d, M_multi=m_multi, min_fidelity=min_fidelity)
    return hp


def predicate_count(dataset):
    info = ROOT / "datasets" / dataset / f"{dataset}.info"
    return sum(1 for l in info.read_text().splitlines() if l.strip())


def completed():
    if not RESULTS.exists():
        return set()
    with open(RESULTS) as f:
        return {(r["dataset"], r["variant"], float(r["epsilon"]), int(r["repeat"])) for r in csv.DictReader(f)}


def run_one(dataset, variant, eps, hp, timeout):
    M = 20 * predicate_count(dataset)
    with tempfile.TemporaryDirectory() as tmp:
        config = json.loads((ROOT / "configs" / f"{dataset}.json").read_text())
        config["output_folder"] = tmp + "/"
        config["force_cyclic"] = variant == "SPECTRUM"
        (Path(tmp) / "psl").mkdir()
        config_path = Path(tmp) / "config.json"
        config_path.write_text(json.dumps(config))
        cmd = [str(BINARIES[variant]), str(config_path), str(M), str(eps), hp["D"], hp["M_multi"], hp["min_fidelity"]]
        env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", str(os.cpu_count())))
        try:
            out = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=timeout, check=True).stdout
        except subprocess.TimeoutExpired:
            return dict(M=M, status="timeout")

    def grab(pattern):
        m = re.search(pattern, out, re.M)
        return m.group(1) if m else ""

    return dict(M=M, random_walks=grab(r"Max number of random walks: (\d+)"),
                fragments=grab(r"^(\d+) Fragments found"), rules=grab(r"Total number of rules (\d+)"),
                runtime_s=grab(r"^SPECTRUM \(Wall Time\) : ([\d.e+-]+) s"), status="ok")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", nargs="+", default=DATASETS)
    parser.add_argument("--epsilons", nargs="+", type=float, default=EPSILONS)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--variants", nargs="+", default=list(BINARIES), choices=list(BINARIES))
    parser.add_argument("--timeout", type=float, default=7200, help="seconds per run")
    args = parser.parse_args()

    for b in BINARIES.values():
        if not b.exists():
            raise SystemExit(f"Missing {b}; run ./build.sh first.")
    hp = hyperparameters()
    done = completed()
    new_file = not RESULTS.exists()
    with open(RESULTS, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            writer.writeheader()
        # cheap runs first: largest epsilon first, all datasets per epsilon
        for eps in sorted(args.epsilons, reverse=True):
            for dataset in args.datasets:
                for variant in args.variants:
                    for repeat in range(1, args.repeats + 1):
                        if (dataset, variant, eps, repeat) in done:
                            continue
                        t0 = time.time()
                        row = dict(dataset=dataset, variant=variant, epsilon=eps, repeat=repeat,
                                   **run_one(dataset, variant, eps, hp[dataset], args.timeout))
                        writer.writerow(row)
                        f.flush()
                        print(f"{dataset:9s} {variant:8s} eps={eps:<6} rep={repeat} fragments={row.get('fragments', '')} "
                              f"rules={row.get('rules', '')} runtime={row.get('runtime_s', '')}s "
                              f"[{row['status']}, {time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
