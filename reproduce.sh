#!/usr/bin/env bash
# Mines rules with SPECTRUM and evaluates them on the test split (NCRL protocol)
# for each dataset in datasets/. Run ./build.sh first.
#
# usage: ./reproduce.sh [family] [umls] [kinship] [fb15k237] [wn18rr]   (default: all five)
#
# Output per dataset in results/<dataset>/:
#   psl/train.json   mined rules
#   spectrum.log     mining log (incl. timings)
#   eval.log         evaluation log (MRR, Hit@1..10)
# and a summary table in results/summary.txt.

set -e
cd "$(dirname "$0")"

SPECTRUM=./build/SPECTRUM
EVALUATOR=./build/spectrum_eval
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(nproc)}

# Hyperparameters per dataset: configs/hyperparameters.tsv; M = 20 x number of predicates in the .info file
declare -A HP=()
while IFS=$'\t' read -r name eps D M_multi min_fidelity; do
  [[ -z "$name" || "$name" == \#* || "$name" == dataset ]] && continue
  HP[$name]="$eps $D $M_multi $min_fidelity"
done < configs/hyperparameters.tsv
M_PER_PREDICATE=20

for bin in "$SPECTRUM" "$EVALUATOR"; do
  [[ -x "$bin" ]] || { echo "Missing $bin; run ./build.sh first." >&2; exit 1; }
done

mkdir -p results
summary=results/summary.txt
printf "%-10s %8s %8s %8s %8s %8s %10s %10s\n" \
  dataset rules MRR Hit@1 Hit@3 Hit@10 mine_s eval_s > "$summary"

for ds in ${@:-family umls kinship fb15k237 wn18rr}; do
  [[ -n "${HP[$ds]}" ]] || { echo "Unknown dataset: $ds" >&2; exit 1; }
  read -r eps D M_multi min_fidelity <<< "${HP[$ds]}"
  M=$(( M_PER_PREDICATE * $(grep -c . "datasets/$ds/$ds.info") ))
  out=results/$ds
  rm -rf "$out"; mkdir -p "$out/psl"

  echo "== $ds: M=$M epsilon=$eps D=$D M_multi=$M_multi min_fidelity=$min_fidelity"
  t0=$(date +%s.%N)
  "$SPECTRUM" "configs/$ds.json" "$M" "$eps" "$D" "$M_multi" "$min_fidelity" > "$out/spectrum.log" 2>&1
  t1=$(date +%s.%N)
  "$EVALUATOR" "datasets/$ds" "$out/psl/train.json" 200 > "$out/eval.log" 2>&1
  t2=$(date +%s.%N)

  read -r mrr h1 h3 h10 < <(awk '/^MRR:/{m=$2} /^Hit@1:/{a=$2} /^Hit@3:/{b=$2} /^Hit@10:/{d=$2} END{print m, a, b, d}' "$out/eval.log")
  rules=$(awk '/^Rules evaluated:/{print $3}' "$out/eval.log")
  printf "%-10s %8s %8s %8s %8s %8s %10.2f %10.2f\n" "$ds" "$rules" "$mrr" "$h1" "$h3" "$h10" \
    "$(echo "$t1 $t0" | awk '{print $1-$2}')" "$(echo "$t2 $t1" | awk '{print $1-$2}')" >> "$summary"
  tail -n 1 "$summary"
done

echo; cat "$summary"
