#!/usr/bin/env bash
# Portable SLURM job for short-read simulations
# ------------------------------------------------------------
#SBATCH --job-name=sim_short
# Keep outputs local to the submission directory (portable default)
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Reasonable defaults (can be overridden at submit time)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=10-00:00:00
# Optional / cluster-specific:
# #SBATCH --partition=general   # <- adjust or comment out as needed

set -euo pipefail
umask 0022

# ======================= CONFIG (override via env) =======================
# Repository root (auto-detected from this file's parent directory)
: "${REPO_ROOT:=$(cd -- "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)}"

# Top-level training root used by this pipeline inside the repo
: "${ROOT:=${REPO_ROOT}/training}"

# Where to put logs produced by the *pipeline* (SLURM stdout/stderr already set above)
: "${PIPELINE_LOG_DIR:=${ROOT}/logs}"

# Tools: discover from PATH, but allow explicit override
: "${ISS_BIN:=$(command -v iss || true)}"
: "${PY:=python3}"

# Simulation parameters (defaults match original script)
: "${MODEL:=miseq}"
: "${SEED:=42}"
: "${CHUNKS:=10}"

# Totals per community
: "${TOTAL_PAIRS_GENERIC:=30000000}"
: "${TOTAL_PAIRS_FILTERED:=8000000}"

# Class ratios (generic and filtered)
: "${GEN_PROK:=0.56}"; : "${GEN_EUK:=0.24}"; : "${GEN_VIR:=0.10}"; : "${GEN_PLA:=0.10}"
: "${FIL_PROK:=0.14}"; : "${FIL_EUK:=0.06}"; : "${FIL_VIR:=0.40}"; : "${FIL_PLA:=0.40}"

# Optional: set SHUFFLE=yes to also write *_shuf_*.fastq.gz
: "${SHUFFLE:=no}"

# Required tools sanity checks
if [[ -z "${ISS_BIN}" || ! -x "${ISS_BIN}" ]]; then
  echo "[ERR] Could not find 'iss'. Set ISS_BIN to the executable path or add to PATH." >&2
  exit 1
fi
if ! command -v md5sum >/dev/null 2>&1; then
  echo "[ERR] md5sum not found in PATH." >&2; exit 1
fi
if ! command -v gzip >/dev/null 2>&1; then
  echo "[ERR] gzip not found in PATH." >&2; exit 1
fi

mkdir -p "${PIPELINE_LOG_DIR}" || true

# ======================= helper functions =======================
timestamp(){ date +%Y%m%d_%H%M%S; }

pick_mate() {
  local prefix="$1" mate="$2"
  for cand in \
    "${prefix}_R${mate}.fastq.gz" "${prefix}_${mate}.fastq.gz" \
    "${prefix}_R${mate}.fastq"    "${prefix}_${mate}.fastq"
  do
    [[ -s "$cand" ]] && { echo "$cand"; return 0; }
  done
  return 1
}

npairs(){ echo $(( $(zcat "$1" | wc -l) / 4 )); }

assert_pairs () {
  local fq1="$1" fq2="$2" want="$3" label="$4"
  local got1 got2
  got1=$(npairs "$fq1"); got2=$(npairs "$fq2")
  [[ "$got1" -eq "$want" && "$got2" -eq "$want" ]] \
    || { echo "[ERR] ${label}: wanted $want pairs, got R1=$got1 R2=$got2"; exit 5; }
}

manifest_json () {
  local label="$1" r1="$2" r2="$3"
  local p1 p2 m1 m2
  p1=$(npairs "$r1"); p2=$(npairs "$r2")
  m1=$(md5sum "$r1" | awk '{print $1}')
  m2=$(md5sum "$r2" | awk '{print $1}')
  cat <<EOF
{"scenario":"$label",
 "r1":{"path":"$r1","pairs":$p1,"md5":"$m1"},
 "r2":{"path":"$r2","pairs":$p2,"md5":"$m2"},
 "model":"$MODEL","seed_base":$SEED,"chunks":$CHUNKS,
 "class_ratios":{
   "prokaryote": $([ "$label" = "GENERIC" ] && echo ${GEN_PROK} || echo ${FIL_PROK}),
   "microeukaryote": $([ "$label" = "GENERIC" ] && echo ${GEN_EUK} || echo ${FIL_EUK}),
   "virus": $([ "$label" = "GENERIC" ] && echo ${GEN_VIR} || echo ${FIL_VIR}),
   "plasmid": $([ "$label" = "GENERIC" ] && echo ${GEN_PLA} || echo ${FIL_PLA})
 }}
EOF
}

shuffle_pairs () {
  local IN1="$1" IN2="$2" OUT1="$3" OUT2="$4"
  paste <(zcat "$IN1" | paste - - - -) <(zcat "$IN2" | paste - - - -) \
    | shuf --random-source=<(yes 42 | head -c 8192) \
    | tee >(cut -f1-4 | tr '\t' '\n' | gzip -c > "$OUT1") \
    | cut -f5-8 | tr '\t' '\n' | gzip -c > "$OUT2"
}

# Per-community helpers defined to use per-community BASE/REFS/ABUND/PERREC/OUT/TMP
concat_class_fasta () {
  local CLASS_LC="$1" CLASS_CAP="$2" REFS="$3" BASE="$4"
  local SEL="${BASE}/genomes/${CLASS_LC}.sel"
  local OUTFA="${REFS}/${CLASS_CAP}.fasta"
  rm -f "$OUTFA"
  [[ -s "$SEL" ]] || { echo "[ERR] Missing selection list: $SEL" >&2; exit 2; }
  while IFS= read -r p; do
    [[ -s "$p" ]] || { echo "[WARN] empty/missing $p"; continue; }
    fname="$(basename "$p")"
    if [[ "$p" == *.gz ]]; then zcat "$p"; else cat "$p"; fi | \
    awk -v grp="$CLASS_CAP" -v fn="$fname" '
      BEGIN{OFS=""}
      /^>/{
        split(substr($0,2), a, /[ \t]/);
        print ">", grp, "|", fn, "|", a[1];
        next
      }
      { print }
    ' >> "$OUTFA"
  done < "$SEL"
  [[ -s "$OUTFA" ]] || { echo "[ERR] Built empty FASTA: $OUTFA" >&2; exit 2; }
}

make_perrec () {
  local CLASS_CAP="$1" REFS="$2" GEN_TSV="$3" OUT_TSV="$4"
  "$PY" - "$REFS/$CLASS_CAP.fasta" "$GEN_TSV" "$OUT_TSV" <<'PY'
import sys, os
from collections import defaultdict
refs, pergen, outp = sys.argv[1], sys.argv[2], sys.argv[3]
g_weight = {}
with open(pergen) as f:
    for line in f:
        if not line.strip(): continue
        parts = line.split()
        if len(parts)<2: continue
        path, w = parts[0], float(parts[1])
        g_weight[os.path.basename(path)] = g_weight.get(os.path.basename(path), 0.0) + w
lens = defaultdict(list)
with open(refs) as f:
    hdr=None; fname=None; L=0
    for line in f:
        if line.startswith('>'):
            if hdr is not None and fname is not None:
                lens[fname].append((hdr, L))
            hdr = line[1:].strip().split()[0]
            parts = hdr.split('|')
            fname = parts[1] if len(parts)>=2 else os.path.basename(hdr)
            L=0
        else:
            L += len(line.strip())
    if hdr is not None and fname is not None:
        lens[fname].append((hdr, L))
with open(outp,'w') as out:
    for fname, recs in lens.items():
        gw = g_weight.get(fname, 0.0)
        tot = sum(L for _,L in recs) or 1.0
        for hdr,L in recs:
            out.write(f"{hdr}\t{gw*(L/tot):.10g}\n")
PY
  [[ -s "$OUT_TSV" ]] || { echo "[ERR] Failed per-record TSV: $OUT_TSV" >&2; exit 2; }
}

gen_group_chunked () {
  local GROUP="$1" REFS_FASTA="$2" ABUND_REC="$3" N_PAIRS="$4" PREFBASE="$5" TMPD="$6"
  local FINAL_R1="${PREFBASE}_${GROUP}_R1.fastq.gz"
  local FINAL_R2="${PREFBASE}_${GROUP}_R2.fastq.gz"
  rm -f "$FINAL_R1" "$FINAL_R2"
  local PER=$(( N_PAIRS / CHUNKS ))
  local REM=$(( N_PAIRS - PER*CHUNKS ))
  echo "[ISS:${GROUP}] pairs=${N_PAIRS} chunks=${CHUNKS} per=${PER} (+1 for first ${REM})" >&2
  for ((i=0;i<CHUNKS;i++)); do
    local Ni_pairs=$PER; if [[ $i -lt $REM ]]; then Ni_pairs=$((Ni_pairs+1)); fi
    [[ $Ni_pairs -gt 0 ]] || continue
    local SEED_I=$(( SEED + i ))
    local PREF="${TMPD}/${GROUP}.c${i}"
    rm -f "${PREF}"_R{1,2}.fastq{,.gz} "${PREF}"_{1,2}.fastq{,.gz} || true
    "$ISS_BIN" generate \
      --genomes "$REFS_FASTA" \
      --abundance_file "$ABUND_REC" \
      --model "$MODEL" \
      --seed "$SEED_I" \
      --n_reads "$((Ni_pairs*2))" \
      --cpus 1 \
      --output "$PREF" \
      --compress
    local R1="$(pick_mate "$PREF" 1 || true)"
    local R2="$(pick_mate "$PREF" 2 || true)"
    [[ -n "$R1" && -n "$R2" ]] || { echo "[ERR] Missing ISS chunk ($GROUP, chunk $i)"; ls -l "${PREF}"* || true; exit 3; }
    if [[ "$R1" == *.gz ]]; then cat "$R1" >> "$FINAL_R1"; else gzip -c "$R1" >> "$FINAL_R1"; fi
    if [[ "$R2" == *.gz ]]; then cat "$R2" >> "$FINAL_R2"; else gzip -c "$R2" >> "$FINAL_R2"; fi
    rm -f "${PREF}"_R{1,2}.fastq{,.gz} "${PREF}"_{1,2}.fastq{,.gz} || true
  done
  gzip -t "$FINAL_R1"; gzip -t "$FINAL_R2"
  echo "$FINAL_R1 $FINAL_R2"
}

process_one() {
  local COMM="$1"
  echo "====== [BEGIN] $COMM ======"
  local BASE="${ROOT}/train_communities/${COMM}"
  local REFS="${BASE}/refs"
  local ABUND="${BASE}/abundance"
  local PERREC="${BASE}/per_record"
  local OUT="${BASE}/reads/short"

  # prefer node-local scratch if available
  local SCRATCH="${SLURM_TMPDIR:-${OUT}}"
  local TMP="${SCRATCH}/.tmp_${SLURM_JOB_ID:-py}_${COMM}"

  mkdir -p "${PIPELINE_LOG_DIR}" "${OUT}" "${TMP}" "${REFS}" "${PERREC}"
  export TMPDIR="${TMP}"; export TMP="${TMP}"

  # archive any pre-existing outputs to avoid accidental mix
  if compgen -G "${OUT}/*.fastq.gz" >/dev/null; then
    local ARCH="${OUT}/.old_$(timestamp)"
    mkdir -p "$ARCH"
    mv "${OUT}"/*.fastq.gz "$ARCH"/ || true
    mv "${OUT}"/*manifest.json "$ARCH"/ 2>/dev/null || true
    echo "[INFO] Archived existing FASTQs to: $ARCH"
  fi

  # Build class FASTAs
  concat_class_fasta prokaryote     Prokaryote     "$REFS" "$BASE"
  concat_class_fasta microeukaryote Microeukaryote "$REFS" "$BASE"
  concat_class_fasta virus          Virus          "$REFS" "$BASE"
  concat_class_fasta plasmid        Plasmid        "$REFS" "$BASE"

  # Per-record weights (from generic TSVs)
  make_perrec Prokaryote     "$REFS" "${ABUND}/prokaryote_generic.tsv"     "${PERREC}/prokaryote_generic_perrecord.tsv"
  make_perrec Microeukaryote "$REFS" "${ABUND}/microeukaryote_generic.tsv" "${PERREC}/microeukaryote_generic_perrecord.tsv"
  make_perrec Virus          "$REFS" "${ABUND}/virus_generic.tsv"          "${PERREC}/virus_generic_perrecord.tsv"
  make_perrec Plasmid        "$REFS" "${ABUND}/plasmid_generic.tsv"        "${PERREC}/plasmid_generic_perrecord.tsv"

  # Class pair counts (use env-configured ratios)
  local N_PROK=$(python3 - <<PY
import math
print(int(round(${TOTAL_PAIRS_GENERIC} * ${GEN_PROK})))
PY
)
  local N_EUK=$(python3 - <<PY
import math
print(int(round(${TOTAL_PAIRS_GENERIC} * ${GEN_EUK})))
PY
)
  local N_VIR=$(python3 - <<PY
import math
print(int(round(${TOTAL_PAIRS_GENERIC} * ${GEN_VIR})))
PY
)
  local N_PLA=$(python3 - <<PY
import math
print(int(round(${TOTAL_PAIRS_GENERIC} * ${GEN_PLA})))
PY
)

  # Generate per-group (GENERIC)
  read -r PROK_R1 PROK_R2 < <(gen_group_chunked "Prokaryote"     "${REFS}/Prokaryote.fasta"     "${PERREC}/prokaryote_generic_perrecord.tsv"     "$N_PROK" "${TMP}/prok" "$TMP")
  read -r EUK_R1  EUK_R2  < <(gen_group_chunked "Microeukaryote" "${REFS}/Microeukaryote.fasta" "${PERREC}/microeukaryote_generic_perrecord.tsv" "$N_EUK"  "${TMP}/euk"  "$TMP")
  read -r VIR_R1  VIR_R2  < <(gen_group_chunked "Virus"          "${REFS}/Virus.fasta"          "${PERREC}/virus_generic_perrecord.tsv"          "$N_VIR"  "${TMP}/vir"  "$TMP")
  read -r PLA_R1  PLA_R2  < <(gen_group_chunked "Plasmid"        "${REFS}/Plasmid.fasta"        "${PERREC}/plasmid_generic_perrecord.tsv"        "$N_PLA"  "${TMP}/pla"  "$TMP")

  local GEN_R1="${OUT}/generic_R1.fastq.gz"
  local GEN_R2="${OUT}/generic_R2.fastq.gz"
  rm -f "$GEN_R1" "$GEN_R2"
  cat "$PROK_R1" "$EUK_R1" "$VIR_R1" "$PLA_R1" > "$GEN_R1"
  cat "$PROK_R2" "$EUK_R2" "$VIR_R2" "$PLA_R2" > "$GEN_R2"
  gzip -t "$GEN_R1"; gzip -t "$GEN_R2"
  assert_pairs "$GEN_R1" "$GEN_R2" "$TOTAL_PAIRS_GENERIC" "GENERIC"
  manifest_json "GENERIC" "$GEN_R1" "$GEN_R2" > "${OUT}/generic_manifest.json"
  echo "[GENERIC] ${COMM}  R1 pairs=$(npairs "$GEN_R1")  R2 pairs=$(npairs "$GEN_R2")"

  if [[ "${SHUFFLE}" == "yes" ]]; then
    shuffle_pairs "$GEN_R1" "$GEN_R2" "${OUT}/generic_shuf_R1.fastq.gz" "${OUT}/generic_shuf_R2.fastq.gz"
    echo "[GENERIC] shuffled copies written"
  fi

  # FILTERED scenario counts
  local F_PROK=$(python3 - <<PY
import math
print(int(round(${TOTAL_PAIRS_FILTERED} * ${FIL_PROK})))
PY
)
  local F_EUK=$(python3 - <<PY
import math
print(int(round(${TOTAL_PAIRS_FILTERED} * ${FIL_EUK})))
PY
)
  local F_VIR=$(python3 - <<PY
import math
print(int(round(${TOTAL_PAIRS_FILTERED} * ${FIL_VIR})))
PY
)
  local F_PLA=$(python3 - <<PY
import math
print(int(round(${TOTAL_PAIRS_FILTERED} * ${FIL_PLA})))
PY
)

  read -r FPROK_R1 FPROK_R2 < <(gen_group_chunked "Prokaryote"     "${REFS}/Prokaryote.fasta"     "${PERREC}/prokaryote_generic_perrecord.tsv"     "$F_PROK" "${TMP}/f_prok" "$TMP")
  read -r FEUK_R1  FEUK_R2  < <(gen_group_chunked "Microeukaryote" "${REFS}/Microeukaryote.fasta" "${PERREC}/microeukaryote_generic_perrecord.tsv" "$F_EUK"  "${TMP}/f_euk" "$TMP")
  read -r FVIR_R1  FVIR_R2  < <(gen_group_chunked "Virus"          "${REFS}/Virus.fasta"          "${PERREC}/virus_generic_perrecord.tsv"          "$F_VIR"  "${TMP}/f_vir" "$TMP")
  read -r FPLA_R1  FPLA_R2  < <(gen_group_chunked "Plasmid"        "${REFS}/Plasmid.fasta"        "${PERREC}/plasmid_generic_perrecord.tsv"        "$F_PLA"  "${TMP}/f_pla" "$TMP")

  local FIL_R1="${OUT}/filtered_R1.fastq.gz"
  local FIL_R2="${OUT}/filtered_R2.fastq.gz"
  rm -f "$FIL_R1" "$FIL_R2"
  cat "$FPROK_R1" "$FEUK_R1" "$FVIR_R1" "$FPLA_R1" > "$FIL_R1"
  cat "$FPROK_R2" "$FEUK_R2" "$FVIR_R2" "$FPLA_R2" > "$FIL_R2"
  gzip -t "$FIL_R1"; gzip -t "$FIL_R2"
  assert_pairs "$FIL_R1" "$FIL_R2" "$TOTAL_PAIRS_FILTERED" "FILTERED"
  manifest_json "FILTERED" "$FIL_R1" "$FIL_R2" > "${OUT}/filtered_manifest.json"
  echo "[FILTERED] ${COMM}  R1 pairs=$(npairs "$FIL_R1")  R2 pairs=$(npairs "$FIL_R2")"

  if [[ "${SHUFFLE}" == "yes" ]]; then
    shuffle_pairs "$FIL_R1" "$FIL_R2" "${OUT}/filtered_shuf_R1.fastq.gz" "${OUT}/filtered_shuf_R2.fastq.gz"
    echo "[FILTERED] shuffled copies written"
  fi

  echo "[OK] ${COMM} short-read simulations complete."
  rm -rf "$TMP"
  echo "====== [END] $COMM ======"
}

# ======================= main =======================
ARG="${1:-ALL}"

if [[ "$ARG" == "ALL" ]]; then
  # loop over all TrainC_short_* communities
  shopt -s nullglob
  for COMM_DIR in "${ROOT}"/train_communities/TrainC_short_*; do
    [[ -d "$COMM_DIR" ]] || continue
    COMM="$(basename "$COMM_DIR")"
    process_one "$COMM"
  done
else
  # single community
  process_one "$ARG"
fi

# ======================= USAGE NOTES =======================
# Submit with defaults:
#   sbatch slurm/sim_short_train_portable.slurm
#
# Run single community:
#   sbatch slurm/sim_short_train_portable.slurm TrainC_short_01
#
# Override parameters at submit time (examples):
#   sbatch --cpus-per-task=4 --mem=64G --time=5-00:00:00 \
#          --export=ALL,TOTAL_PAIRS_GENERIC=2000000,CHUNKS=4,SHUFFLE=yes \
#          slurm/sim_short_train_portable.slurm
#
# Use a conda env or modules before submission, or add to the top:
#   module load python/3.10
#   module load iss/2.0
#   # or
#   source ~/miniconda3/etc/profile.d/conda.sh && conda activate myenv
