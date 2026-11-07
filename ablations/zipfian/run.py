#!/usr/bin/env python3
"""Zipfian / top-M assumption: are useful but low-frequency rules missed by SPECTRUM's mining?

For each dataset this mines twice, keeping every candidate rule (M_multi = --keep, so get_final_rules
truncates nothing) with the same min_fidelity, so only the mining budget differs:

  reference: N (paths per node) = --reference-N, by solving N = M ln M / (eps^2 |V|) for eps with the
             heuristic M. With N far above the number of paths in any node's D-neighbourhood, mining is
             exhaustive, so this is the set of rules SPECTRUM can learn at all.
  default:   the dataset's hyperparameters (configs/hyperparameters.tsv), --repeats times (mining samples
             paths at random wherever N is smaller than a neighbourhood).
  reference_all: the reference run with min_fidelity 0, giving the exhaustive precision of every pattern.
             Used to check rules a default run admits but the reference rejects ("spurious"): their
             sampled precision is above min_fidelity while their exhaustive precision is not.

Each reference rule's precision ("fidelity": precision x symmetry factor / Bayesian prior), recall
(sum over the head facts it covers of ln(1 + groundings entailing that fact)) and utility come from
psl/train_stats.csv. Utility includes precision, and is normalised there to sum to 1 over the rules with
the same head predicate, so rarity is measured by recall. The script records, per reference rule, its recall
and precision percentiles, in how
many default runs it was mined, and in how many it was also in the returned rule set (the top M_multi
candidates by precision for its head predicate, which is what SPECTRUM outputs).

min_fidelity defaults to the dataset's value in configs/hyperparameters.tsv (0.1 for family and wn18rr,
1 = precision better than chance for the others), as in reproduce.sh. It filters candidate rules after mining, so it
does not change what is mined.

Rules are matched across runs by pattern: head predicate and multiset of body predicates. A pattern
does not record edge directions, and print_rule takes the variable bindings from the first stored
location, so where sampling stores a different location first the same pattern prints with different
variable order (and, being sampled, a different precision); each pattern prints as exactly one rule. The
stricter exact match (variable names and body order canonicalised) is kept in the *_exact columns.

Outputs: ablations/zipfian/rules.csv (one row per reference rule), summary.csv (one row per dataset);
raw rule files go to ablations/zipfian/out/<dataset>/<run>/ and are reused on a rerun.

usage: python3 ablations/zipfian/run.py [--datasets ...] [--repeats 3] [--reference-N 1e9] [--min-fidelity F]
Then: python3 ablations/zipfian/plot.py
"""
import argparse
import csv
import itertools
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import ROOT, check_binaries, heuristic_M, hyperparameters, run_miner  # noqa: E402

OUT = HERE / "out"
DATASETS = ["family", "umls", "kinship", "wn18rr", "fb15k237"]
INT_MAX = 2**31 - 1  # main.cpp stores N in an int


def number_of_nodes(dataset):
    """Distinct constants in train.db (HyperGraph::number_of_nodes)."""
    nodes = set()
    for line in (ROOT / "datasets" / dataset / "train.db").read_text().splitlines():
        if "(" in line and not line.startswith("//"):
            nodes.update(a.strip() for a in line[line.index("(") + 1:line.rindex(")")].split(","))
    return len(nodes)


def atoms(text):
    """[(predicate, [args])] for 'p(A,B) & q(B,C)'. Predicate names may contain '&' (UMLS); print_rule
    writes a predicate's trailing digits as a quoted constant argument, e.g. term(A,B,'21') (Kinship)."""
    out = []
    for atom in text.strip().split(" & "):
        name, args = atom[:atom.index("(")], atom[atom.index("(") + 1:atom.rindex(")")]
        out.append((name, [a.strip() for a in args.split(",")]))
    return out


def canonical(rule):
    """Rule string up to variable names and body-atom order: 'b1(A,B) & b2(B,C) >> h(A,C)'.
    Quoted constants are kept as they are."""
    body, head = rule.split(" >> ")
    head_atom, body_atoms = atoms(head)[0], atoms(body)
    best = None
    for perm in itertools.permutations(body_atoms):
        names = {}

        def write(atom):
            name, args = atom
            return name + "(" + ",".join(a if a.startswith("'") else names.setdefault(a, f"V{len(names)}")
                                         for a in args) + ")"
        h = write(head_atom)
        s = " & ".join(write(a) for a in perm) + " >> " + h
        best = s if best is None or s < best else best
    return best


def head_predicate(rule):
    """Head predicate, including a constant split off by print_rule: 'term' + '17' -> term'17'."""
    name, args = atoms(rule.split(" >> ")[1])[0]
    return name + "".join(a for a in args if a.startswith("'"))


def pattern(rule):
    """(head predicate, sorted body predicates): the rule's pattern, ignoring variable bindings."""
    body, head = rule.split(" >> ")
    key = lambda a: a[0] + "".join(x for x in a[1] if x.startswith("'"))  # noqa: E731
    return key(atoms(head)[0]), tuple(sorted(key(a) for a in atoms(body)))


def read_stats(folder):
    with open(folder / "psl" / "train_stats.csv") as f:
        return [dict(rule=r["rule"], utility=float(r["utility"]), fidelity=float(r["fidelity"]),
                     recall=float(r["recall"]))
                for r in csv.DictReader(f)]


def mine(dataset, name, M, eps, hp, timeout):
    """Mines once into out/<dataset>/<name>/ unless that run is already there; returns its run info."""
    folder = OUT / dataset / name
    info_file = folder / "run.json"
    if info_file.exists():
        return json.loads(info_file.read_text()), folder
    res = run_miner("SPECTRUM", dataset, M, f"{eps:.10g}", hp, timeout, details=True, output_folder=folder)
    info = dict(M=M, epsilon=eps, **res)
    if res["status"] == "ok":
        info_file.write_text(json.dumps(info, indent=1))
    print(f"{dataset:9s} {name:9s} N={info.get('random_walks', '')} candidates={info.get('rules', '')} "
          f"runtime={info.get('runtime_s', '')}s [{info['status']}]", flush=True)
    return info, folder


def missed_rate(rows, lo, hi, repeats):
    """Mean percentage of the rules with recall percentile in (lo, hi] missed per default run."""
    sel = [r for r in rows if lo < float(r["recall_pct"]) <= hi]
    return 100.0 * sum(repeats - r["mined_runs"] for r in sel) / (repeats * len(sel)) if sel else 0.0


def percentile(values, v):
    """Percentage of values <= v."""
    return 100.0 * sum(x <= v for x in values) / len(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", nargs="+", default=DATASETS)
    parser.add_argument("--repeats", type=int, default=3, help="default-hyperparameter runs per dataset")
    parser.add_argument("--reference-N", type=float, default=1e9, help="paths per node for the reference run")
    parser.add_argument("--min-fidelity", default=None,
                        help="candidate precision threshold for both runs (default: the dataset's own value)")
    parser.add_argument("--keep", default="1000000", help="M_multi for both runs (keep every candidate)")
    parser.add_argument("--timeout", type=float, default=21600, help="seconds per run")
    args = parser.parse_args()
    if args.reference_N > INT_MAX:
        raise SystemExit(f"--reference-N must be at most {INT_MAX} (N is an int in main.cpp)")

    check_binaries(["SPECTRUM"])
    hp_all = hyperparameters()
    rule_rows, summary_rows = [], []
    for dataset in args.datasets:
        min_fidelity = args.min_fidelity if args.min_fidelity is not None else hp_all[dataset]["min_fidelity"]
        hp = dict(hp_all[dataset], M_multi=args.keep, min_fidelity=min_fidelity)
        M, V = heuristic_M(dataset), number_of_nodes(dataset)
        eps_ref = math.sqrt(M * math.log(M) / (args.reference_N * V))
        ref_info, ref_folder = mine(dataset, "reference", M, eps_ref, hp, args.timeout)
        all_info, all_folder = mine(dataset, "reference_all", M, eps_ref, dict(hp, min_fidelity="0"), args.timeout)
        defaults = [mine(dataset, f"default_{i}", M, hp_all[dataset]["epsilon"], hp, args.timeout)
                    for i in range(1, args.repeats + 1)]
        if ref_info["status"] != "ok" or all_info["status"] != "ok" or any(info["status"] != "ok" for info, _ in defaults):
            print(f"{dataset}: a run did not finish; skipping its analysis", flush=True)
            continue

        reference = read_stats(ref_folder)
        recalls = [r["recall"] for r in reference]
        fidelities = [r["fidelity"] for r in reference]
        m_multi = int(hp_all[dataset]["M_multi"])
        mined, returned, mined_exact = defaultdict(int), defaultdict(int), defaultdict(int)
        default_patterns = set()
        exhaustive = {pattern(r["rule"]): r["fidelity"] for r in read_stats(all_folder)}
        reference_patterns = {pattern(r["rule"]) for r in reference}
        spurious = []  # exhaustive precision of default rules the reference rejects
        for _, folder in defaults:
            rules = read_stats(folder)
            by_head = defaultdict(list)
            for r in rules:
                mined[pattern(r["rule"])] += 1
                if pattern(r["rule"]) not in reference_patterns:
                    spurious.append(exhaustive.get(pattern(r["rule"]), float("nan")))
                mined_exact[canonical(r["rule"])] += 1
                default_patterns.add(pattern(r["rule"]))
                by_head[head_predicate(r["rule"])].append(r)
            for head_rules in by_head.values():  # what SPECTRUM returns: top M_multi by precision per head
                for r in sorted(head_rules, key=lambda r: -r["fidelity"])[:m_multi]:
                    returned[pattern(r["rule"])] += 1

        for r in reference:
            p = pattern(r["rule"])
            rule_rows.append(dict(dataset=dataset, rule=r["rule"], recall=r["recall"], fidelity=r["fidelity"],
                                  utility=r["utility"], recall_pct=f"{percentile(recalls, r['recall']):.1f}",
                                  fidelity_pct=f"{percentile(fidelities, r['fidelity']):.1f}",
                                  mined_runs=mined[p], returned_runs=returned[p],
                                  mined_runs_exact=mined_exact[canonical(r["rule"])], repeats=args.repeats))
        n = len(reference)
        rows = rule_rows[-n:]
        summary_rows.append(dict(
            dataset=dataset, reference_N=ref_info["random_walks"], default_N=defaults[0][0]["random_walks"],
            min_fidelity=min_fidelity, reference_rules=n,
            mined_every_run=sum(r["mined_runs"] == args.repeats for r in rows),
            missed_some_run=sum(r["mined_runs"] < args.repeats for r in rows),
            missed_every_run=sum(r["mined_runs"] == 0 for r in rows),
            # mean fraction missed per default run, by recall quartile (q1 = lowest recall)
            **{f"missed_pct_recall_q{q + 1}": f"{missed_rate(rows, 25 * q, 25 * (q + 1), args.repeats):.1f}"
               for q in range(4)},
            returned_every_run=sum(r["returned_runs"] == args.repeats for r in rows),
            missed_some_run_exact=sum(r["mined_runs_exact"] < args.repeats for r in rows),
            spurious_per_run=f"{len(spurious) / args.repeats:.1f}",
            spurious_exhaustive_precision_median=f"{sorted(spurious)[len(spurious) // 2]:.3f}" if spurious else "",
            spurious_exhaustive_precision_max=f"{max(spurious):.3f}" if spurious else ""))
        s = summary_rows[-1]
        print(f"{dataset}: {n} reference rules; mined in every default run: {s['mined_every_run']}, "
              f"missed in some run: {s['missed_some_run']}, returned every run: {s['returned_every_run']}", flush=True)

    done = {r["dataset"] for r in summary_rows}
    for name, rows in (("rules.csv", rule_rows), ("summary.csv", summary_rows)):
        if (HERE / name).exists():  # keep the rows of datasets not analysed in this invocation
            with open(HERE / name) as f:
                rows = [r for r in csv.DictReader(f) if r["dataset"] not in done] + rows
        if rows:
            with open(HERE / name, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
    print(f"Wrote {HERE / 'rules.csv'} and {HERE / 'summary.csv'}")


if __name__ == "__main__":
    main()
