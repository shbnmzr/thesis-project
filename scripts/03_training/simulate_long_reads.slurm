#!/usr/bin/env bash
# Portable SLURM job for long-read simulations with NanoSim
# ------------------------------------------------------------
#SBATCH --job-name=sim_long_train
# Keep outputs local to the submission directory
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Parallelism & resources (override at submit time if needed)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=20-00:00:00
# Optional / cluster-specific:
# #SBATCH --partition=general

set -euo pipefail
umask 0022

# ======================= INPUT ARG =======================
COMM="${1:-}"
if [[ -z "$COMM" ]]; then
  echo "Usage: sbatch slurm/simulate_long_reads.slurm <TrainC_long_XX>" >&2
  exit 1
fi

# ======================= CONFIG (override via env) =======================
# Repository root (auto-detected from this file's parent directory)
: "${REPO_ROOT:=$(cd -- "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)}"

# Top-level training root used by this pipeline inside the repo
: "${ROOT:=${REPO_ROOT}/training}"

# Expose BASE so the embedded Python can see it
export BASE="${ROOT}/train_communities/${COMM}"

# Data & output
ABUND="${BASE}/abundance"
OUT="${BASE}/reads/long"
LOGDIR="${ROOT}/logs"

# Temporary workspace (local scratch if available, else project dir)
SCRATCH="${SLURM_TMPDIR:-${OUT}}"
TMP="${SCRATCH}/.tmp_${SLURM_JOB_ID:-py}_${COMM}"

mkdir -p "$OUT" "$TMP" "$LOGDIR"

# ======================= TOOLS =======================
# Prefer PATH discovery; allow explicit override via env
: "${PYTHON_BIN:=$(command -v python3 || true)}"
: "${SEQTK:=$(command -v seqtk || true)}"
: "${PARALLEL_BIN:=$(command -v parallel || true)}"

# NanoSim installation (either SIMULATOR_PY path or `nanosim` module)
# If you use a source checkout, set SIMULATOR_PY to path/to/src/simulator.py
: "${SIMULATOR_PY:=}"

# Sanity checks
[[ -n "${PYTHON_BIN}" && -x "${PYTHON_BIN}" ]] || { echo "[ERR] python3 not found (set PYTHON_BIN)." >&2; exit 2; }
[[ -n "${SEQTK}" && -x "${SEQTK}" ]] || { echo "[ERR] seqtk not found (set SEQTK or add to PATH)." >&2; exit 2; }
[[ -n "${PARALLEL_BIN}" && -x "${PARALLEL_BIN}" ]] || { echo "[ERR] GNU parallel not found (set PARALLEL_BIN or add to PATH)." >&2; exit 2; }

if [[ -z "${SIMULATOR_PY}" ]]; then
  # Try to locate NanoSim simulator.py via Python import path
  SIMULATOR_PY="$(${PYTHON_BIN} - <<'PY'
import sys, importlib.util
spec = importlib.util.find_spec('nanosim')
if not spec or not spec.submodule_search_locations:
    print('')
else:
    import os
    # common layout .../nanosim/src/simulator.py (pip wheels may differ)
    for root in spec.submodule_search_locations:
        cand = os.path.join(root, 'src', 'simulator.py')
        if os.path.exists(cand):
            print(cand)
            break
PY
)"
fi
[[ -n "${SIMULATOR_PY}" && -s "${SIMULATOR_PY}" ]] || { echo "[ERR] NanoSim simulator.py not found. Set SIMULATOR_PY to the file path." >&2; exit 2; }

# ======================= NanoSim model =======================
# Option A: provide a tarball to unpack once per community
: "${NANOSIM_MODEL_TAR:=}"
# Option B: point to an existing prepared model directory containing a 'model' subdir
: "${MODEL_DIR:=${BASE}/.nanosim_model}"

if [[ ! -d "${MODEL_DIR}/model" ]]; then
  if [[ -n "${NANOSIM_MODEL_TAR}" && -s "${NANOSIM_MODEL_TAR}" ]]; then
    mkdir -p "${MODEL_DIR}/model"
    tar -xzf "${NANOSIM_MODEL_TAR}" -C "${MODEL_DIR}/model" --strip-components=1
  else
    echo "[ERR] NanoSim model not found. Either set NANOSIM_MODEL_TAR to a valid tar.gz or provide MODEL_DIR with a 'model' subdirectory." >&2
    exit 2
  fi
fi

# ======================= Targets & ratios =======================
# Long-read totals per scenario (override if desired)
: "${TOTAL_READS_GENERIC:=300000}"
: "${TOTAL_READS_FILTERED:=80000}"

# Class ratios (generic and filtered); defaults mirror your originals
: "${GEN_PROK_RATIO:=0.56}"; : "${GEN_EUK_RATIO:=0.24}"; : "${GEN_VIR_RATIO:=0.10}"; : "${GEN_PLA_RATIO:=0.10}"
: "${FIL_PROK_RATIO:=0.14}"; : "${FIL_EUK_RATIO:=0.06}"; : "${FIL_VIR_RATIO:=0.40}"; : "${FIL_PLA_RATIO:=0.40}"

# Compute integer targets per class
calc_target(){ ${PYTHON_BIN} - "$@" <<'PY'
import sys, math
T = int(sys.argv[1])
ratios = list(map(float, sys.argv[2:]))
# largest remainder to hit exact total
floors = [int(math.floor(T*r)) for r in ratios]
rem = T - sum(floors)
order = sorted(range(len(ratios)), key=lambda i: (T*ratios[i]-floors[i]), reverse=True)
for i in order[:rem]:
    floors[i]+=1
print('\n'.join(map(str,floors)))
PY
}
read GEN_PROK GEN_EUK GEN_VIR GEN_PLA < <(calc_target "$TOTAL_READS_GENERIC" "$GEN_PROK_RATIO" "$GEN_EUK_RATIO" "$GEN_VIR_RATIO" "$GEN_PLA_RATIO")
read FIL_PROK FIL_EUK FIL_VIR FIL_PLA < <(calc_target "$TOTAL_READS_FILTERED" "$FIL_PROK_RATIO" "$FIL_EUK_RATIO" "$FIL_VIR_RATIO" "$FIL_PLA_RATIO")

# We'll simulate ONCE per class at MAX(Generic, Filtered) to reuse pools
export SIM_PROK=$(( GEN_PROK > FIL_PROK ? GEN_PROK : FIL_PROK ))
export SIM_EUK=$(( GEN_EUK > FIL_EUK ? GEN_EUK : FIL_EUK ))
export SIM_VIR=$(( GEN_VIR > FIL_VIR ? GEN_VIR : FIL_VIR ))
export SIM_PLA=$(( GEN_PLA > FIL_PLA ? GEN_PLA : FIL_PLA ))

# Prune tiny per-genome allocations and redistribute (speeds up dramatically)
: "${MIN_READS_PER_GENOME:=50}"; export MIN_READS_PER_GENOME

echo "[PLAN] ${COMM}: simulate class pools (long reads)"
printf "  Prokaryote     pool: %d\n" "$SIM_PROK"
printf "  Microeukaryote pool: %d\n" "$SIM_EUK"
printf "  Virus          pool: %d\n" "$SIM_VIR"
printf "  Plasmid        pool: %d\n" "$SIM_PLA"

# ======================= build per-class plan =======================
PLAN="${TMP}/plan.tsv"
/usr/bin/env python3 - "$PLAN" <<'PY'
import sys, math, os
from pathlib import Path

base = Path(os.environ["BASE"])
out_plan = Path(sys.argv[1])

SIM_PROK = int(os.environ["SIM_PROK"])
SIM_EUK  = int(os.environ["SIM_EUK"])
SIM_VIR  = int(os.environ["SIM_VIR"])
SIM_PLA  = int(os.environ["SIM_PLA"])
MINR     = int(os.environ.get("MIN_READS_PER_GENOME", "50"))

tsv = {
 "Prokaryote":     base/"abundance/prokaryote_generic.tsv",
 "Microeukaryote": base/"abundance/microeukaryote_generic.tsv",
 "Virus":          base/"abundance/virus_generic.tsv",
 "Plasmid":        base/"abundance/plasmid_generic.tsv",
}

target = {"Prokaryote":SIM_PROK,"Microeukaryote":SIM_EUK,"Virus":SIM_VIR,"Plasmid":SIM_PLA}

def load_items(tsvfile):
    items=[]
    with open(tsvfile) as f:
        for ln in f:
            if not ln.strip(): continue
            p,w = ln.split()[:2]
            items.append((p, float(w)))
    return items

with out_plan.open('w') as out:
    for G in ("Prokaryote","Microeukaryote","Virus","Plasmid"):
        items = load_items(tsv[G])
        T = target[G]
        if T <= 0 or not items: continue

        # Initial largest-remainder rounding
        s = sum(w for _,w in items) or 1.0
        floors=[]; rema=[]; sf=0
        for p,w in items:
            r = (w/s) * T
            k = int(math.floor(r))
            floors.append([p,k]); rema.append((r-k, p)); sf += k
        rem = T - sf
        for _,p in sorted(rema, reverse=True)[:rem]:
            for row in floors:
                if row[0]==p: row[1]+=1; break

        # Prune tiny allocations and redistribute proportionally
        small = sum(k for _,k in floors if k < MINR)
        keep = [(p,k) for p,k in floors if k >= MINR]
        if small > 0 and keep:
            s_keep = sum(k for _,k in keep)
            adj = []
            carry = small
            for p,k in keep:
                add = int((k / s_keep) * small)
                adj.append([p, k + add]); carry -= add
            for i in sorted(range(len(adj)), key=lambda i: keep[i][1], reverse=True)[:max(carry,0)]:
                adj[i][1] += 1
            floors = adj
        else:
            floors = keep if keep else floors

        for p,k in floors:
            if k>0: out.write(f"{p}\t{k}\t{G}\n")
PY

echo "[PLAN] $(wc -l < "$PLAN") genome rows"

# ======================= simulate per genome in parallel =======================
mkdir -p "${TMP}/Prokaryote" "${TMP}/Microeukaryote" "${TMP}/Virus" "${TMP}/Plasmid"

export PYTHON_BIN SIMULATOR_PY MODEL_DIR TMP
"${PARALLEL_BIN}" -j "${SLURM_CPUS_PER_TASK}" --colsep '\t' '
  i={#}; GENOME={1}; READS={2}; GRP={3}
  [[ "$READS" -gt 0 ]] || exit 0
  PREF="${TMP}/${GRP}/ns_${i}"
  (
    cd "${MODEL_DIR}/model" && \
    "${PYTHON_BIN}" "${SIMULATOR_PY}" genome -r "$GENOME" -n "$READS" -o "$PREF" --seed $((42+i))
  ) && (
    [[ -s "${PREF}_reads.fastq" ]] && mv "${PREF}_reads.fastq" "${TMP}/${GRP}/chunk_${i}.fastq" \
      || [[ -s "${PREF}_reads.fasta" ]] && mv "${PREF}_reads.fasta" "${TMP}/${GRP}/chunk_${i}.fastq" \
      || echo "[WARN] No reads for ${GENOME}"
  )
' :::: "$PLAN"

# ======================= Merge class pools =======================
POOL_PROK="${OUT}/pool_prokaryote.fastq"
POOL_EUK="${OUT}/pool_microeukaryote.fastq"
POOL_VIR="${OUT}/pool_virus.fastq"
POOL_PLA="${OUT}/pool_plasmid.fastq"

if compgen -G "${TMP}/Prokaryote/chunk_*.fastq" > /dev/null; then cat "${TMP}/Prokaryote"/chunk_*.fastq > "$POOL_PROK"; else : > "$POOL_PROK"; fi
if compgen -G "${TMP}/Microeukaryote/chunk_*.fastq" > /dev/null; then cat "${TMP}/Microeukaryote"/chunk_*.fastq > "$POOL_EUK"; else : > "$POOL_EUK"; fi
if compgen -G "${TMP}/Virus/chunk_*.fastq" > /dev/null; then cat "${TMP}/Virus"/chunk_*.fastq > "$POOL_VIR"; else : > "$POOL_VIR"; fi
if compgen -G "${TMP}/Plasmid/chunk_*.fastq" > /dev/null; then cat "${TMP}/Plasmid"/chunk_*.fastq > "$POOL_PLA"; else : > "$POOL_PLA"; fi

# ======================= helpers =======================
nreads(){ echo $(( $(wc -l < "$1") / 4  )); }

sample_exact () {
  local IN="$1" TARGET="$2" OUT="$3" LABEL="$4"
  "${SEQTK}" sample -s 42 "$IN" "$TARGET" > "$OUT"
  local GOT=$(( $(wc -l < "$OUT") / 4 ))
  if [[ "$GOT" -ne "$TARGET" ]]; then
    echo "[ERR] ${LABEL}: wanted $TARGET reads, got $GOT" >&2
    exit 4
  fi
}

# Sanity log
echo "[POOL] counts:"
printf "  prok: %s\n" "$(nreads "$POOL_PROK")"
printf "  euk : %s\n" "$(nreads "$POOL_EUK")"
printf "  vir : %s\n" "$(nreads "$POOL_VIR")"
printf "  pla : %s\n" "$(nreads "$POOL_PLA")"

# ======================= derive Generic & Filtered =======================
GENERIC="${OUT}/generic.fastq"
FILTERED="${OUT}/filtered.fastq"
rm -f "$GENERIC" "$FILTERED"

# Generic targets (parallel)
GEN_PROK_F="${TMP}/gen_prok.fastq"; sample_exact "$POOL_PROK" "$GEN_PROK" "$GEN_PROK_F" "prok(generic)" &
GEN_EUK_F="${TMP}/gen_euk.fastq";  sample_exact "$POOL_EUK"  "$GEN_EUK"  "$GEN_EUK_F"  "euk(generic)"  &
GEN_VIR_F="${TMP}/gen_vir.fastq";  sample_exact "$POOL_VIR"  "$GEN_VIR"  "$GEN_VIR_F"  "vir(generic)"  &
GEN_PLA_F="${TMP}/gen_pla.fastq";  sample_exact "$POOL_PLA"  "$GEN_PLA"  "$GEN_PLA_F"  "pla(generic)"  &
wait
cat "$GEN_PROK_F" "$GEN_EUK_F" "$GEN_VIR_F" "$GEN_PLA_F" > "$GENERIC"

# Filtered targets (parallel)
FIL_PROK_F="${TMP}/fil_prok.fastq"; sample_exact "$POOL_PROK" "$FIL_PROK" "$FIL_PROK_F" "prok(filtered)" &
FIL_EUK_F="${TMP}/fil_euk.fastq";  sample_exact "$POOL_EUK"  "$FIL_EUK"  "$FIL_EUK_F"  "euk(filtered)"  &
FIL_VIR_F="${TMP}/fil_vir.fastq";  sample_exact "$POOL_VIR"  "$FIL_VIR"  "$FIL_VIR_F"  "vir(filtered)"  &
FIL_PLA_F="${TMP}/fil_pla.fastq";  sample_exact "$POOL_PLA"  "$FIL_PLA"  "$FIL_PLA_F"  "pla(filtered)"  &
wait
cat "$FIL_PROK_F" "$FIL_EUK_F" "$FIL_VIR_F" "$FIL_PLA_F" > "$FILTERED"

echo "[DONE] ${COMM} long-read simulations"
echo "  Generic : ${GENERIC}  (reads=$(nreads "$GENERIC"))"
echo "  Filtered: ${FILTERED} (reads=$(nreads "$FILTERED"))"

# Optional: compress final outputs to save space (uncomment if desired)
# if command -v pigz >/dev/null 2>&1; then pigz -p ${SLURM_CPUS_PER_TASK} "$GENERIC" "$FILTERED"; else gzip "$GENERIC" "$FILTERED"; fi

# Cleanup big temps
rm -rf "$TMP"

# ======================= USAGE NOTES =======================
# Submit with defaults:
#   sbatch slurm/simulate_long_reads.slurm TrainC_long_01
#
# Override resources at submit time:
#   sbatch --cpus-per-task=32 --mem=64G --time=10-00:00:00 \
#          slurm/simulate_long_reads.slurm TrainC_long_01
#
# Override model location and totals via env export:
#   sbatch --export=ALL,NANOSIM_MODEL_TAR=/path/to/model.tar.gz,TOTAL_READS_GENERIC=200000 \
#          slurm/simulate_long_reads.slurm TrainC_long_01
#
# Use a conda env or modules before submission, or add to the top:
#   module load python/3.10
#   module load nanosim/3
#   module load seqtk
#   # or
#   source ~/miniconda3/etc/profile.d/conda.sh && conda activate nanosim
