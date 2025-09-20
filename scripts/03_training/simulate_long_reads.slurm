#!/bin/bash
#SBATCH --job-name=sim_long_train_%x
#SBATCH --output=/path/to/training/logs/%x_%A.out
#SBATCH --error=/path/to/training/logs/%x_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=20-00:00:00
#SBATCH --partition=allgroups

set -euo pipefail
umask 0022

# Locale & threads
export LC_ALL=C
export LANG=C
export LC_CTYPE=C
export OMP_NUM_THREADS=1

COMM="${1:?Need community id, e.g., TrainC_long_01}"

# ---- paths ----
ROOT="/path/to/training"
export BASE="${ROOT}/train_communities/${COMM}"
ABUND="${BASE}/abundance"
OUT="${BASE}/reads/long"
TMP="${OUT}/.tmp_${SLURM_JOB_ID}"
LOGDIR="${ROOT}/logs"
mkdir -p "$OUT" "$TMP" "$LOGDIR"

# ---- tools (NanoSim env only) ----
PYTHON="/path/to/miniconda3/envs/nanosim/bin/python"
SEQTK="/path/to//miniconda3/envs/nanosim/bin/seqtk"
NANOSIM_REPO="/path/to/tools/NanoSim"
SIMULATOR_PY="${NANOSIM_REPO}/src/simulator.py"

command -v parallel >/dev/null 2>&1 || { echo "[ERR] GNU parallel not found"; exit 2; }
[[ -x "$PYTHON" ]] || { echo "[ERR] Python not found at $PYTHON"; exit 2; }
[[ -x "$SEQTK"  ]] || { echo "[ERR] seqtk not found at $SEQTK"; exit 2; }
[[ -s "$SIMULATOR_PY" ]] || { echo "[ERR] NanoSim simulator.py missing at $SIMULATOR_PY"; exit 2; }

# ---- Check scikit-learn version ----
REQ_SKLEARN="0.23.2"
CUR_SKLEARN="$("$PYTHON" - <<'PY'
import sklearn, sys
print(sklearn.__version__)
PY
)"
if [[ "$CUR_SKLEARN" != "$REQ_SKLEARN" ]]; then
  echo "[ERR] scikit-learn version is ${CUR_SKLEARN}, but model requires ${REQ_SKLEARN}. Please pin your nanosim env to scikit-learn==${REQ_SKLEARN}." >&2
  exit 9
fi

# ---- NanoSim model (cached per community) ----
MODEL_TAR="/nfsd/bcb/bcbg/tools/NanoSim/pre-trained_models/metagenome_ERR3152364_Even_v3.2.2.tar.gz"
MODEL_DIR="${BASE}/.nanosim_model"
if [[ ! -d "${MODEL_DIR}/model" ]]; then
  mkdir -p "${MODEL_DIR}/model"
  tar -xzf "$MODEL_TAR" -C "${MODEL_DIR}/model" --strip-components=1
fi

MODEL_PREFIX="training"
export MODEL_PREFIX
echo "[INFO] NanoSim model prefix: $MODEL_PREFIX"

# ---- targets (per community) ----
# Generic total: 300k reads => prok 168k, euk 72k, virus 30k, plasmid 30k
# Filtered total: 80k reads => prok 11.2k, euk 4.8k, virus 32k, plasmid 32k
GEN_PROK=168000; GEN_EUK=72000; GEN_VIR=30000; GEN_PLA=30000
FIL_PROK=11200;  FIL_EUK=4800;  FIL_VIR=32000;  FIL_PLA=32000

# Simulate ONCE per class at MAX(Generic, Filtered).
export SIM_PROK=$GEN_PROK
export SIM_EUK=$GEN_EUK
export SIM_VIR=$FIL_VIR
export SIM_PLA=$FIL_PLA

# Prune tiny per-genome allocations and redistribute
MIN_READS_PER_GENOME="${MIN_READS_PER_GENOME:-500}"
export MIN_READS_PER_GENOME

echo "[PLAN] ${COMM}: simulate class pools (long reads)"
printf "  Prokaryote     pool: %d\n" "$SIM_PROK"
printf "  Microeukaryote pool: %d\n" "$SIM_EUK"
printf "  Virus          pool: %d\n" "$SIM_VIR"
printf "  Plasmid        pool: %d\n" "$SIM_PLA"

# ---- build per-class plans (largest-remainder rounding + prune) ----
PLAN="${TMP}/plan.tsv"
"$PYTHON" - "$PLAN" <<'PY'
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

# ---- simulate per genome in parallel; produce per-class POOLS ----
mkdir -p "${TMP}/Prokaryote" "${TMP}/Microeukaryote" "${TMP}/Virus" "${TMP}/Plasmid"

PAR_OPTS="--jobs ${SLURM_CPUS_PER_TASK} --line-buffer --delay 0.2 --halt now,fail=1"
export PYTHON SIMULATOR_PY MODEL_DIR TMP MODEL_PREFIX
parallel ${PAR_OPTS} --colsep '\t' '
  i={#}; GENOME={1}; READS={2}; GRP={3}
  [[ "$READS" -gt 0 ]] || exit 0
  PREF="${TMP}/${GRP}/ns_${i}"
  (
    cd "$MODEL_DIR/model" \
    && "$PYTHON" "$SIMULATOR_PY" genome -r "$GENOME" -n "$READS" -o "$PREF" \
       --seed $((42+i)) -t 1 --fastq -min 50 -c "$MODEL_PREFIX"
  )
  rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[ERR] NanoSim failed (rc=$rc) for $GENOME" >&2
    exit $rc
  fi
  if [[ -s "${PREF}_reads.fastq" ]]; then
    mv "${PREF}_reads.fastq" "${TMP}/${GRP}/chunk_${i}.fastq"
  elif [[ -s "${PREF}_reads.fasta" ]]; then
    mv "${PREF}_reads.fasta" "${TMP}/${GRP}/chunk_${i}.fastq"
  else
    echo "[WARN] No reads for ${GENOME}"
  fi
' :::: "$PLAN"

# ---- Merge class pools ----
POOL_PROK="${OUT}/pool_prokaryote.fastq"
POOL_EUK="${OUT}/pool_microeukaryote.fastq"
POOL_VIR="${OUT}/pool_virus.fastq"
POOL_PLA="${OUT}/pool_plasmid.fastq"

if compgen -G "${TMP}/Prokaryote/chunk_*.fastq" > /dev/null; then cat "${TMP}/Prokaryote"/chunk_*.fastq > "$POOL_PROK"; else : > "$POOL_PROK"; fi
if compgen -G "${TMP}/Microeukaryote/chunk_*.fastq" > /dev/null; then cat "${TMP}/Microeukaryote"/chunk_*.fastq > "$POOL_EUK"; else : > "$POOL_EUK"; fi
if compgen -G "${TMP}/Virus/chunk_*.fastq" > /dev/null; then cat "${TMP}/Virus"/chunk_*.fastq > "$POOL_VIR"; else : > "$POOL_VIR"; fi
if compgen -G "${TMP}/Plasmid/chunk_*.fastq" > /dev/null; then cat "${TMP}/Plasmid"/chunk_*.fastq > "$POOL_PLA"; else : > "$POOL_PLA"; fi

# ---- Require pools to be large enough before sampling ----
need_pool () {
  local file="$1" label="$2" need="$3"
  local have=0
  [[ -s "$file" ]] && have=$(( $(wc -l < "$file") / 4 ))
  if [[ ! -s "$file" ]]; then
    echo "[ERR] ${label}: pool file missing/empty: $file" >&2
    exit 5
  fi
  if (( have < need )); then
    echo "[ERR] ${label}: pool has ${have} reads, need ${need}. Upstream NanoSim produced too few reads." >&2
    exit 6
  fi
}

need_pool "$POOL_PROK" "prok pool" "$GEN_PROK"
need_pool "$POOL_EUK"  "euk pool"  "$GEN_EUK"
need_pool "$POOL_VIR"  "virus pool" "$GEN_VIR"
need_pool "$POOL_PLA"  "plasmid pool" "$GEN_PLA"

# ---- helpers ----
nreads(){ echo $(( $(wc -l < "$1") / 4  )); }

sample_exact () {
  local IN="$1" TARGET="$2" OUT="$3" LABEL="$4"
  "$SEQTK" sample -s 42 "$IN" "$TARGET" > "$OUT"
  local GOT=$(( $(wc -l < "$OUT") / 4 ))
  if [[ "$GOT" -ne "$TARGET" ]]; then
    echo "[ERR] ${LABEL}: wanted $TARGET reads, got $GOT" >&2
    exit 4
  fi
}

# ---- derive Generic & Filtered from pools (parallel sampling) ----
GENERIC="${OUT}/generic.fastq"
FILTERED="${OUT}/filtered.fastq"
rm -f "$GENERIC" "$FILTERED"

GEN_PROK_F="${TMP}/gen_prok.fastq"; sample_exact "$POOL_PROK" "$GEN_PROK" "$GEN_PROK_F" "prok(generic)" &
GEN_EUK_F="${TMP}/gen_euk.fastq";  sample_exact "$POOL_EUK"  "$GEN_EUK"  "$GEN_EUK_F"  "euk(generic)"  &
GEN_VIR_F="${TMP}/gen_vir.fastq";  sample_exact "$POOL_VIR"  "$GEN_VIR"  "$GEN_VIR_F"  "vir(generic)"  &
GEN_PLA_F="${TMP}/gen_pla.fastq";  sample_exact "$POOL_PLA"  "$GEN_PLA"  "$GEN_PLA_F"  "pla(generic)"  &
wait
cat "$GEN_PROK_F" "$GEN_EUK_F" "$GEN_VIR_F" "$GEN_PLA_F" > "$GENERIC"

FIL_PROK_F="${TMP}/fil_prok.fastq"; sample_exact "$POOL_PROK" "$FIL_PROK" "$FIL_PROK_F" "prok(filtered)" &
FIL_EUK_F="${TMP}/fil_euk.fastq";  sample_exact "$POOL_EUK"  "$FIL_EUK"  "$FIL_EUK_F"  "euk(filtered)"  &
FIL_VIR_F="${TMP}/fil_vir.fastq";  sample_exact "$POOL_VIR"  "$FIL_VIR"  "$FIL_VIR_F"  "vir(filtered)"  &
FIL_PLA_F="${TMP}/fil_pla.fastq";  sample_exact "$POOL_PLA"  "$FIL_PLA"  "$FIL_PLA_F"  "pla(filtered)"  &
wait
cat "$FIL_PROK_F" "$FIL_EUK_F" "$FIL_VIR_F" "$FIL_PLA_F" > "$FILTERED"

echo "[POOL] counts:"
printf "  prok: %s\n" "$(nreads "$POOL_PROK")"
printf "  euk : %s\n" "$(nreads "$POOL_EUK")"
printf "  vir : %s\n" "$(nreads "$POOL_VIR")"
printf "  pla : %s\n" "$(nreads "$POOL_PLA")"

echo "[DONE] ${COMM} long-read simulations"
echo "  Generic : ${GENERIC}  (reads=$(nreads "$GENERIC"))"
echo "  Filtered: ${FILTERED} (reads=$(nreads "$FILTERED"))"

# Compress final outputs to save space
pigz -p ${SLURM_CPUS_PER_TASK} "$GENERIC" "$FILTERED"

# Cleanup temps
rm -rf "$TMP"
