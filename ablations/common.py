"""Shared helpers for the ablation studies: hyperparameters, running a miner, results CSVs."""
import csv
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARIES = {"SPECTRUM": ROOT / "build" / "SPECTRUM", "RW": ROOT / "build" / "SPECTRUM_RW"}


def hyperparameters():
    """Per-dataset hyperparameters from configs/hyperparameters.tsv."""
    hp = {}
    for line in (ROOT / "configs" / "hyperparameters.tsv").read_text().splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("dataset"):
            continue
        name, eps, d, m_multi, min_fidelity = line.split("\t")
        hp[name] = dict(epsilon=float(eps), D=d, M_multi=m_multi, min_fidelity=min_fidelity)
    return hp


def heuristic_M(dataset):
    """M = 20 x number of predicates in the dataset's .info file."""
    info = ROOT / "datasets" / dataset / f"{dataset}.info"
    return 20 * sum(1 for l in info.read_text().splitlines() if l.strip())


def check_binaries(variants):
    for v in variants:
        if not BINARIES[v].exists():
            raise SystemExit(f"Missing {BINARIES[v]}; run ./build.sh first.")


def run_miner(variant, dataset, M, eps, hp, timeout, data_folder=None, details=False, output_folder=None):
    """Runs SPECTRUM (pattern mining) or SPECTRUM_RW (random walks, force_cyclic off) once.
    Returns random walks, fragments found, candidate rules mined and mining wall time.
    data_folder replaces the dataset's train.db folder (same .info file); details=True also returns
    the pattern locations generated before deduplication and the fragment-creation (mining) time.
    output_folder keeps the rule files (psl/train.json, psl/train_stats.csv) there instead of discarding them."""
    with tempfile.TemporaryDirectory() as tmp:
        if output_folder is not None:
            Path(output_folder).mkdir(parents=True, exist_ok=True)
            tmp = str(output_folder)
        config = json.loads((ROOT / "configs" / f"{dataset}.json").read_text())
        config["output_folder"] = tmp + "/"
        if data_folder is not None:
            config["data_folder"] = str(data_folder).rstrip("/") + "/"
        config["force_cyclic"] = variant == "SPECTRUM"
        (Path(tmp) / "psl").mkdir(exist_ok=True)
        config_path = Path(tmp) / "config.json"
        config_path.write_text(json.dumps(config))
        cmd = [str(BINARIES[variant]), str(config_path), str(M), str(eps), hp["D"], hp["M_multi"], hp["min_fidelity"]]
        env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", str(os.cpu_count())))
        try:
            out = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=timeout, check=True).stdout
        except subprocess.TimeoutExpired:
            return dict(status="timeout")

    def grab(pattern):
        m = re.search(pattern, out, re.M)
        return m.group(1) if m else ""

    res = dict(random_walks=grab(r"Max number of random walks: (\d+)"),
               fragments=grab(r"^(\d+) Fragments found"), rules=grab(r"Total number of rules (\d+)"),
               runtime_s=grab(r"^SPECTRUM \(Wall Time\) : ([\d.e+-]+) s"), status="ok")
    if details:
        res.update(locations=grab(r"^Before deduplication location count:(\d+)"),
                   mining_s=grab(r"^Fragment Creation \(Wall Time\) : ([\d.e+-]+) s"))
    return res


def plain_log_ticks(ax):
    """Plain-number tick labels on a log y axis (190, 200, 1000 rather than 1.9x10^2). Call after
    plotting. Axes spanning less than a factor of 10 get evenly spaced labelled ticks, since a log
    axis would otherwise label at most one or two values there."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, LogFormatter, MaxNLocator

    class Plain(LogFormatter):
        def _num_to_string(self, x, vmin, vmax):
            return f"{x:g}"

    lo, hi = ax.get_ylim()
    ax.yaxis.set_major_formatter(Plain(labelOnlyBase=False, minor_thresholds=(2, 0.4)))
    if hi / lo < 10:
        ticks = [t for t in MaxNLocator(5).tick_values(lo, hi) if lo <= t <= hi and t > 0]
        ax.yaxis.set_minor_locator(FixedLocator(ticks))
        ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    else:
        ax.yaxis.set_minor_formatter(Plain(labelOnlyBase=False, minor_thresholds=(2, 0.4)))


def read_results(path):
    if not Path(path).exists():
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


class ResultsWriter:
    """Appends rows to a results CSV, writing the header for a new file."""

    def __init__(self, path, fields):
        new = not Path(path).exists()
        self.f = open(path, "a", newline="")
        self.writer = csv.DictWriter(self.f, fieldnames=fields)
        if new:
            self.writer.writeheader()

    def write(self, row):
        self.writer.writerow(row)
        self.f.flush()

    def close(self):
        self.f.close()
