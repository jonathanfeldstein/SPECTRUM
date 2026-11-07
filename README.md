# SPECTRUM

SPECTRUM is an optimized C++ rule mining algorithm that learns logical rules from knowledge graphs and relational databases. 

It works by first mining pattern fragments in the data and then using these fragments to compute a heuristic utility of candidate rules. Finally, candidate rules are then optimised and ranked so as to maximise the utility of the set of rules.

## Table of Contents
- [Quick start: reproducing the results](#quick-start-reproducing-the-results)
- [Repository layout](#repository-layout)
- [Dependencies](#dependencies)
- [Building the Project](#building-the-project)
- [Running SPECTRUM](#running-spectrum)
- [Evaluating rules](#evaluating-rules)
- [Code Overview](#code-overview)
- [Reproducing Ablation Studies](#reproducing-ablation-studies)

## Quick start: reproducing the results

```bash
./build.sh
./reproduce.sh                 # all five datasets, or e.g. ./reproduce.sh family umls
```

Afterwards, `python3 plot_hits.py` plots Hit@1–Hit@10 for every dataset (`results/hits_at_k.png` / `.pdf`; needs matplotlib).

For each dataset this mines rules from `datasets/<dataset>/train.db` and evaluates them on the
test split. Rules and logs go to `results/<dataset>/`, and a summary to `results/summary.txt`.
Expected output:

| Dataset | Rules evaluated | MRR | Hit@1 | Hit@3 | Hit@10 | Mining (s) | Evaluation (s) |
|---|---|---|---|---|---|---|---|
| family | 160 | 0.9011 | 0.8194 | 0.9820 | 0.9993 | 0.3 | 0.2 |
| umls | 11 | 0.7606 | 0.5905 | 0.9143 | 0.9524 | 0.4 | 0.1 |
| kinship | 752 | 0.4514 | 0.2771 | 0.5425 | 0.8335 | 0.4 | 0.3 |
| wn18rr | 41 | 0.5004 | 0.3315 | 0.5916 | 0.8899 | 1.1 | 0.5 |
| fb15k237 | ~2960 | ~0.501 | ~0.406 | ~0.548 | ~0.682 | 18–24 | 6.5 |

family, umls, kinship and wn18rr are deterministic and reproduce these numbers exactly. On fb15k237 the
random walks used for fragment mining are seeded non-deterministically, so the mined rules and the
metrics vary slightly between runs (about ±0.001 in Hit@k).

The hyperparameters are set in `configs/hyperparameters.tsv`. M is 20 × the number of predicates listed in
`datasets/<dataset>/<dataset>.info`:

| Dataset | Predicates | M | ε | D | M_multi | min_fidelity |
|---|---|---|---|---|---|---|
| family | 12 | 240 | 0.01 | 3 | 20 | 0.1 |
| umls | 44 | 880 | 0.01 | 3 | 20 | 1 |
| kinship | 24 | 480 | 0.01 | 3 | 40 | 1 |
| wn18rr | 11 | 220 | 0.01 | 3 | 20 | 0.1 |
| fb15k237 | 237 | 4740 | 0.01 | 3 | 20 | 1 |
| cad | 19 | 380 | 0.01 | 3 | 20 | 1 |
| yelp | 5 | 100 | 0.01 | 3 | 20 | 1 |

SPECTRUM mines from `datasets/<dataset>/train.db` only. For the five KGC datasets the evaluator then scores
`test.txt`, following rule bodies over the background facts, training and validation triples (NCRL's
protocol). CAD and Yelp have no test split: they are only mined (`ablations/scaling`). See
[datasets/README.md](datasets/README.md) for which split is used for what, the sources, and how CAD's
and Yelp's `train.db` were built from the original PSL files (`datasets/prepare_psl_datasets.py`).

## Repository layout

```
build.sh, CMakeLists.txt   build SPECTRUM, SPECTRUM_RW and the evaluator into build/
reproduce.sh               mine + evaluate every dataset → results/
plot_hits.py               Hit@1–10 figure from results/
src/                       SPECTRUM rule miner
evaluator/                 rule evaluator (NCRL protocol)
configs/<dataset>.json     mining configuration per dataset
configs/hyperparameters.tsv  SPECTRUM hyperparameters per dataset
datasets/<dataset>/        mining data (train.db, .info) and evaluation data (.txt); see datasets/README.md
dependencies/              bundled Boost, Eigen, Spectra
ablations/                 ablation studies (see below)
exp/                       older ablation scripts
```

## Dependencies

- **OS**: Linux (macOS with Homebrew is also supported by the build).
- **C++ Compiler**: supporting C++20 (e.g., GCC 10+ or Clang 10+).
- **CMake**: 3.16 or higher.
- **OpenMP**: for multi-threading.
- **jemalloc**: e.g. `apt install libjemalloc-dev`, or install it into an active conda/micromamba
  environment (`build.sh` adds `$CONDA_PREFIX` to CMake's search path).
- **Bundled Libraries** in `dependencies/`, linked automatically:
  - Boost (1.83.0)
  - Eigen (3.4.0)
  - Spectra (1.0.1)
- **Python 3.x**: only for the ablation scripts.

Tested with GCC 11.4 and CMake 3.25 on Ubuntu 22.04.

## Building the Project

Run

```bash
./build.sh
```

This will create a `build/` directory containing:
- `SPECTRUM`: the rule mining executable.
- `SPECTRUM_RW`: SPECTRUM with pattern mining replaced by random walks (used in the ablation study); same command line.
- `spectrum_eval`: the rule evaluator.

## Running SPECTRUM

```bash
./build/SPECTRUM <config_file.json> <M> <epsilon> <D> <M_multi> [min_fidelity=1]
```

- `<config_file.json>`: dataset paths, target predicates, output folder and logic (see `configs/`).
- `<M>`: number of candidate rules to consider; SPECTRUM selects an optimised subset of the `M` rules with highest heuristic utility.
- `<epsilon>`: permitted upper bound on the uncertainty in rule utility (smaller epsilon requires more computation).
- `<D>`: maximum rule length.
- `<M_multi>`: maximum number of rules output per head predicate.
- `[min_fidelity]`: minimum fidelity for a candidate rule to be kept (default 1).

The output folder named in the config must contain a `psl/` subfolder; rules are written to
`<output_folder>/psl/train.json`. For example:

```bash
mkdir -p results/family/psl
./build/SPECTRUM configs/family.json 240 0.01 3 20 0.1
```

## Evaluating rules

```bash
./build/spectrum_eval datasets/<dataset> results/<dataset>/psl/train.json [top_k=200]
```

The evaluator implements the knowledge-graph completion protocol of NCRL (Cheng et al., ICLR 2023)
and gives results identical to NCRL's `kg_completion.py` on the same rules. It is 10–170× faster.
For each test triple `a r b` it ranks all entities as candidates for `a`:

- Only chain-shaped rules (`r(X,Y) ← b1(X,Z1) ∧ … ∧ bn(Zn-1,Y)`, with atoms in either direction) are
  used, at most `top_k` per head relation, in the order SPECTRUM outputs them. Rules whose
  relations do not occur in `train.txt` or whose relation names contain `&` are not used.
- A candidate's score is the sum over rules of (rule weight × number of body paths from `b` to
  the candidate), counted over `facts.txt` + inverse facts + `train.txt` + `valid.txt`.
- Known true answers other than `a` are filtered out of the ranking.
- Test queries whose relation has no rule are not scored; the evaluator prints how many were.

- The rank of `a` is 1 + the number of candidates scoring strictly higher than `a`, as in NCRL.
  The evaluator prints MRR and Hit@1–Hit@10 over the scored queries.

## Code Overview

- **`src/main.cpp`**: The main entry point for SPECTRUM.
- **`src/HyperGraphs/`**: Contains the core data structures (`HyperGraph`, `Relation`) used to represent relational data and datasets in memory.
- **`src/Fragments/`**: Implements the logic for generating, storing, and managing rule fragments (both standard cyclic/acyclic and random-walk specific).
- **`src/CandidateRules/`**: Manages the construction, proposal, and validation of final candidate logical rules mined from the graph.
- **`src/Utils/`**: A collection of utility helpers for file I/O operations (`read_files_utils`), performance profiling (`profiling_utils`), and data structure helpers (`map_utils`, `set_utils`, `vector_utils`).
- **`evaluator/`**: `Graph.hpp` loads a dataset folder, `RuleParser.hpp` turns SPECTRUM's PSL rules into relation paths, `Evaluator.hpp` scores and ranks the test queries.

## Reproducing Ablation Studies

The project contains automated Python scripts to reproduce the ablation studies for `M` and `epsilon` (`ablations/M_epsilon`) and for pattern mining vs random walks (`ablations/random_walk_comparison`). 

**1. SPECTRUM as M and epsilon vary:**
```bash
python3 ablations/M_epsilon/run.py    # appends to ablations/M_epsilon/results.csv
python3 ablations/M_epsilon/plot.py   # → ablations/M_epsilon/figure.png / figure.pdf (needs matplotlib, numpy)
```

M ablation: ε = 0.01 and M = 2, 4, 8, … up to 16 × the heuristic M (20 × number of predicates), plus the
heuristic M itself, marked with a star in the runtime panel. ε ablation: M = heuristic M and
ε ∈ {0.1, 0.05, 0.01, 0.005, 0.001}. 3 repeats per setting; by default on UMLS, Kinship, Family and WN18RR
(`--datasets` adds e.g. `fb15k237`). Like the other ablation, `run.py` resumes from `results.csv`.

**2. Pattern mining vs random walks as a function of epsilon:**
```bash
python3 ablations/random_walk_comparison/run.py  # runs SPECTRUM and SPECTRUM_RW; appends to ablations/random_walk_comparison/results.csv
python3 ablations/random_walk_comparison/plot.py # → ablations/random_walk_comparison/figure.png / figure.pdf (needs matplotlib, numpy)
```

For each dataset, ε ∈ {0.5, 0.1, 0.05, 0.01, 0.005} and 3 repeats, `run.py` records the number of
fragments (patterns) found, candidate rules mined and the mining wall time of both miners. All other
hyperparameters are the dataset's values in `configs/hyperparameters.tsv`. `run.py` skips runs already
in `results.csv`, so it can be stopped and resumed, and `--datasets`, `--epsilons`, `--repeats` select
a subset. The random-walk runs at small ε are slow on FB15K-237 (about 0.5 h per run at ε = 0.01 and
2 h at ε = 0.005). The figure plots rules and runtime only: the random-walk miner runs
with `force_cyclic` off, so its fragment counts include non-cyclic fragments that SPECTRUM does not
count (they remain in `results.csv`).

**3. Runtime vs dataset size by replication:**
```bash
python3 ablations/replication/run.py   # appends to ablations/replication/results.csv
python3 ablations/replication/plot.py  # → ablations/replication/figure.png / figure.pdf (needs matplotlib, numpy)
```

Each dataset is copied k ∈ {1, 2, 4, 8, 16} times as disjoint graphs, so nodes, facts and every local
neighbourhood grow by exactly k. With ε_k = ε/√k the number of paths per node N stays fixed and the
mining analysis predicts runtime ∝ k; with ε fixed (N ∝ 1/k) the total path budget M ln M / ε² stays
constant. Single-threaded by default (`--threads`), 3 repeats, UMLS, Kinship, Family and WN18RR
(about 20 minutes). Panel (c) plots runtime against the pattern locations generated.

**4. Are useful but low-frequency rules missed? (Zipfian / top-M assumption):**
```bash
python3 ablations/zipfian/run.py       # → ablations/zipfian/rules.csv, summary.csv (raw rule files in ablations/zipfian/out/)
python3 ablations/zipfian/plot.py      # → ablations/zipfian/figure.png / figure.pdf (needs matplotlib)
```

For each dataset, a reference run with N = 10⁹ paths per node (exhaustive mining) lists every rule
SPECTRUM can learn with precision × symmetry / prior above 1, with its utility; three runs with the
default hyperparameters show which of them are mined, and which are in the returned rule set (top
`M_multi` per head by precision). Every candidate is kept in both (`M_multi` = 10⁶), so only the
mining budget differs. Rules are matched by pattern (head and body predicates), since
the variable order printed for a pattern depends on which location was stored first. A second reference run without the
precision cut-off checks default-run rules whose sampled precision is above 1 but whose exhaustive
precision is not. The figure plots precision against recall (the rarity of a rule; utility is not used
because it includes precision) and marks rules missed in a default run. The miner writes recall as the
last column of `_stats.csv`. About 5 minutes; FB15K-237's
reference runs take about a minute each.

**5. Mining time vs dataset size:**
```bash
python3 ablations/scaling/run.py       # appends to ablations/scaling/results.csv
python3 ablations/scaling/plot.py      # → ablations/scaling/figure.png / figure.pdf (needs matplotlib, numpy)
```

Mines all seven datasets (including CAD and Yelp) 5 times each with the hyperparameters above and plots the
mining wall time against the number of facts in `train.db`, with a log-log fit. Run it on an
otherwise idle machine. About 5 minutes, most of it FB15K-237, CAD and Yelp.
