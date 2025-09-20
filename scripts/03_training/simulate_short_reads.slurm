#!/usr/bin/env bash
#SBATCH --job-name=sim_short_train
#SBATCH --output=/path/to/training/logs/%x_%A.out
#SBATCH --error=/path/to/training/logs/%x_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=10-00:00:00
#SBATCH --partition=allgroups

set -euo pipefail
umask 0022
ulimit -f unlimited 2>/dev/null || true

# ---- paths & tools ----
ROOT="/path/to/training"
ISS_BIN="/path/to/miniconda3/envs/python39/bin/iss"
PY="/path/to/miniconda3/envs/python39/bin/python3"

# ---- sim params ----
MODEL="${MODEL:-hiseq}"
SEED=${SEED:-42}
CHUNKS=${CHUNKS:-10}

# totals per community
TOTAL_PAIRS_GENERIC=${TOTAL_PAIRS_GENERIC:-30000000}
TOTAL_PAIRS_FILTERED=${TOTAL_PAIRS_FILTERED:-8000000}

# class ratios (documentary only; used in manifest)
GEN_PROK=0.56; GEN_EUK=0.24; GEN_VIR=0.10; GEN_PLA=0.10
FIL_PROK=0.14; FIL_EUK=0.06; FIL_VIR=0.40; FIL_PLA=0.40

# options
SHUFFLE="${SHUFFLE:-no}"
VERIFY_CHUNKS="${VERIFY_CHUNKS:-yes}"
CHUNK_MAX_RETRIES="${CHUNK_MAX_RETRIES:-10}"

timestamp(){ date +%Y%m%d_%H%M%S; }

# --- sanity: tools present ---
need_tool(){ command -v "$1" >/dev/null 2>&1 || { echo "[ERR] Missing required tool: $1" >&2; exit 1; }; }
need_tool gzip; need_tool awk; need_tool paste; need_tool shuf; need_tool md5sum
[[ -x "$ISS_BIN" ]] || { echo "[ERR] iss not found/executable at: $ISS_BIN" >&2; exit 1; }
[[ -x "$PY" ]] || { echo "[ERR] python not found/executable at: $PY" >&2; exit 1; }

npairs(){
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    gzip -dc "$f" | wc -l | awk '{print int($1/4)}'
  else
    wc -l < "$f" | awk '{print int($1/4)}'
  fi
}

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

assert_pairs () {
  local fq1="$1" fq2="$2" want="$3" label="$4"
  local got1 got2
  got1=$(npairs "$fq1"); got2=$(npairs "$fq2")
  if [[ "$got1" -ne "$want" || "$got2" -ne "$want" ]]; then
    echo "[ERR] ${label}: wanted $want pairs, got R1=$got1 R2=$got2" >&2
    exit 5
  fi
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
   "prokaryote": $([ "$label" = "GENERIC" ] && echo 0.56 || echo 0.14),
   "microeukaryote": $([ "$label" = "GENERIC" ] && echo 0.24 || echo 0.06),
   "virus": $([ "$label" = "GENERIC" ] && echo 0.10 || echo 0.40),
   "plasmid": $([ "$label" = "GENERIC" ] && echo 0.10 || echo 0.40)
 }}
EOF
}

shuffle_pairs () {
  local IN1="$1" IN2="$2" OUT1="$3" OUT2="$4"
  paste <(gzip -dc "$IN1" | paste - - - -) <(gzip -dc "$IN2" | paste - - - -) \
    | shuf --random-source=<(yes "$SEED") \
    | tee >(cut -f1-4 | tr '\t' '\n' | gzip -c > "$OUT1") \
    | cut -f5-8 | tr '\t' '\n' | gzip -c > "$OUT2"
}

gz_repack_concat() {
  local out="$1"; shift
  {
    for f in "$@"; do
      if [[ "$f" == *.gz ]]; then gzip -dc "$f"; else cat "$f"; fi
    done
  } | gzip -1c > "${out}.part"
  gzip -t "${out}.part"
  mv -f "${out}.part" "$out"
}

atom_mv() {
  local src="$1" dst="$2"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --no-owner --no-group "$src" "${dst}.part" && mv -f "${dst}.part" "$dst"
  else
    cp -f "$src" "${dst}.part" && mv -f "${dst}.part" "$dst"
  fi
  [[ -s "$dst" ]] || { echo "[ERR] move failed: $src -> $dst" >&2; exit 6; }
}

check_space() {
  local path="$1" need_gb="$2"
  local avail mp
  mp="$(df -P "$path" | awk 'NR==2{print $6}')"
  avail="$(df -P "$path" | awk 'NR==2{print $4}')" # 1K blocks
  local avail_gb=$(( avail / 1024 / 1024 ))
  echo "[INFO] Free space on $mp: ~${avail_gb}G"
  if [[ "$avail_gb" -lt "$need_gb" ]]; then
    echo "[WARN] Low free space on $mp (~${avail_gb}G < ${need_gb}G)." >&2
  fi
}

# ---- Per-community helpers ----
concat_class_fasta () {
  local CLASS_LC="$1" CLASS_CAP="$2" REFS="$3" BASE="$4"
  local SEL="${BASE}/genomes/${CLASS_LC}.sel"
  local OUTFA="${REFS}/${CLASS_CAP}.fasta"
  rm -f "$OUTFA"
  [[ -s "$SEL" ]] || { echo "[ERR] Missing selection list: $SEL" >&2; exit 2; }
  while IFS= read -r p; do
    [[ -s "$p" ]] || { echo "[WARN] empty/missing $p" >&2; continue; }
    fname="$(basename "$p")"
    if [[ "$p" == *.gz ]]; then gzip -dc "$p"; else cat "$p"; fi | \
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

# AWK filter for short contigs (>=500 bp by default; adjust via editing min=)
filter_fasta_minlen () {
  local IN="$1" OUT="$2" MIN="$3"
  awk -v min="$MIN" '
  BEGIN{hdr=""; seq=""}
  /^>/{
    if (hdr != "" && length(seq) >= min) { print hdr; print seq; }
    hdr = $0; seq = ""; next
  }
  { gsub(/\r/,""); seq = seq $0 }
  END{ if (hdr != "" && length(seq) >= min) { print hdr; print seq; } }
  ' "$IN" > "$OUT"
  [[ -s "$OUT" ]] || { echo "[ERR] Filtering produced empty FASTA: $OUT" >&2; exit 7; }
}

make_perrec () {
  local CLASS_CAP="$1" REFS="$2" GEN_TSV="$3" OUT_TSV="$4"
  # normalized per-record weights; writes TSV summing to 1.0
  "$PY" - "$REFS/$CLASS_CAP.fasta" "$GEN_TSV" "$OUT_TSV" <<'PY'
import sys, os
from collections import defaultdict

refs, pergen, outp = sys.argv[1], sys.argv[2], sys.argv[3]

# Read genome-level weights (may sum to < 1 if class is a fraction of a larger mix)
g_weight = {}
with open(pergen) as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        path, w = parts[0], float(parts[1])
        g_weight[os.path.basename(path)] = g_weight.get(os.path.basename(path), 0.0) + w

# Map contigs (records) to their source filename and length
lens = defaultdict(list)
with open(refs) as f:
    hdr=None; fname=None; L=0
    for line in f:
        if line.startswith('>'):
            if hdr is not None and fname is not None:
                lens[fname].append((hdr, L))
            hdr = line[1:].strip().split()[0]
            parts = hdr.split('|')
            fname = parts[1] if len(parts) >= 2 else os.path.basename(hdr)
            L = 0
        else:
            L += len(line.strip())
    if hdr is not None and fname is not None:
        lens[fname].append((hdr, L))

# Distribute genome weights to per-record weights
records = []  # (hdr, raw_weight)
for fname, recs in lens.items():
    gw = g_weight.get(fname, 0.0)
    tot = sum(L for _, L in recs) or 1.0
    for hdr, L in recs:
        records.append((hdr, gw * (L / tot)))

# Normalize per-record weights to sum to 1.0 exactly
S = sum(w for _, w in records)
with open(outp, 'w') as out:
    if S <= 0.0:
        # Fallback: uniform
        n = len(records) or 1
        for hdr, _ in records:
            out.write("{}\t{:.10g}\n".format(hdr, 1.0/n))
    else:
        invS = 1.0 / S
        for hdr, w in records:
            out.write("{}\t{:.10g}\n".format(hdr, w*invS))
PY
  [[ -s "$OUT_TSV" ]] || { echo "[ERR] Failed per-record TSV: $OUT_TSV" >&2; exit 2; }
}

# ---- chunked generator (never uses --paired; tops up until exact) ----
gen_group_chunked () {
  local GROUP="$1" REFS_FASTA="$2" ABUND_REC="$3" N_PAIRS="$4" PREFBASE="$5" TMPD="$6"

  local FINAL_R1="${PREFBASE}_${GROUP}_R1.fastq.gz"
  local FINAL_R2="${PREFBASE}_${GROUP}_R2.fastq.gz"
  rm -f "$FINAL_R1" "$FINAL_R2"

  local PER=$(( N_PAIRS / CHUNKS ))
  local REM=$(( N_PAIRS - PER*CHUNKS ))
  echo "[ISS:${GROUP}] pairs=${N_PAIRS} chunks=${CHUNKS} per=${PER} (+1 for first ${REM})" >&2

  local PIECES_R1=()
  local PIECES_R2=()

  for ((i=0;i<CHUNKS;i++)); do
    local want=$PER; (( i < REM )) && want=$((want+1))
    [[ $want -gt 0 ]] || continue

    local produced=0
    local attempts=0

    while [[ $produced -lt $want ]]; do
      attempts=$((attempts+1))
      if [[ $attempts -gt $CHUNK_MAX_RETRIES ]]; then
        echo "[ERR] ${GROUP} chunk c${i}: could not reach ${want} pairs after ${CHUNK_MAX_RETRIES} attempts (stuck at ${produced})" >&2
        exit 4
      fi
      local need=$(( want - produced ))
      local SEED_I=$(( SEED + i*100 + attempts ))
      local PREF="${TMPD}/${GROUP}.c${i}.t${attempts}"
      rm -f "${PREF}"_R{1,2}.fastq{,.gz} "${PREF}"_{1,2}.fastq{,.gz} || true

      "$ISS_BIN" generate \
        --genomes "$REFS_FASTA" \
        --abundance_file "$ABUND_REC" \
        --model "$MODEL" \
        --seed "$SEED_I" \
        --n_reads "$((need*2))" \
        --cpus 1 \
        --output "$PREF" \
        --compress

      local R1="$(pick_mate "$PREF" 1 || true)"
      local R2="$(pick_mate "$PREF" 2 || true)"
      [[ -n "$R1" && -n "$R2" ]] || { echo "[ERR] Missing ISS output ($GROUP c${i} t${attempts})" >&2; ls -l "${PREF}"* || true; exit 3; }

      local got1 got2
      got1=$(npairs "$R1"); got2=$(npairs "$R2")
      [[ "$got1" -eq "$got2" ]] || { echo "[ERR] Mate mismatch ($GROUP c${i} t${attempts}): R1=$got1 R2=$got2" >&2; exit 4; }

      if [[ "$VERIFY_CHUNKS" == "yes" ]]; then
        local ts; ts="$(date +%F_%T)"
        if [[ -n "${CHK_LOG:-}" ]]; then
          printf "%s\t%s\tc%d\tt%d\tneed=%d\tgot=%d\ttotal=%d\twant=%d\n" \
            "$ts" "$GROUP" "$i" "$attempts" "$need" "$got1" "$((produced+got1))" "$want" >> "$CHK_LOG"
        else
          echo "CHK $GROUP c${i} t${attempts}: need ${need}, got ${got1}, total $((produced+got1)) / ${want}" >&2
        fi
      fi

      produced=$(( produced + got1 ))
      PIECES_R1+=("$R1")
      PIECES_R2+=("$R2")
    done
  done

  gz_repack_concat "$FINAL_R1" "${PIECES_R1[@]}"
  gz_repack_concat "$FINAL_R2" "${PIECES_R2[@]}"
  gzip -t "$FINAL_R1"; gzip -t "$FINAL_R2"

  if [[ "$VERIFY_CHUNKS" == "yes" ]]; then
    local got_final1 got_final2
    got_final1=$(npairs "$FINAL_R1")
    got_final2=$(npairs "$FINAL_R2")
    if [[ $got_final1 -ne $N_PAIRS || $got_final2 -ne $N_PAIRS ]]; then
      echo "[ERR] ${GROUP} final: want $N_PAIRS pairs, got R1=$got_final1 R2=$got_final2" >&2
      exit 4
    fi
  fi

  rm -f "${PIECES_R1[@]}" "${PIECES_R2[@]}"
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

  # temp dir selection
  local TMPROOT="${SLURM_TMPDIR:-}"
  if [[ -z "$TMPROOT" || ! -w "$TMPROOT" || ! -x "$TMPROOT" ]]; then
    TMPROOT="$OUT"
    echo "[INFO] Using NFS for temp (no writable SLURM_TMPDIR)."
  else
    echo "[INFO] Using node-local SLURM_TMPDIR for temp: $TMPROOT"
  fi
  TMP="${TMPROOT}/.tmp_${SLURM_JOB_ID:-$$}_${COMM}"
  mkdir -p "${ROOT}/logs" "$OUT" "$TMP" "$REFS" "$PERREC"
  export TMPDIR="${TMP}"

  # cleanup trap
  cleanup(){ rc=$?; [[ -n "${TMP:-}" && -d "$TMP" ]] && rm -rf "$TMP"; exit $rc; }
  trap cleanup EXIT INT TERM

  check_space "$OUT" 50 || true

  # archive existing outputs
  if compgen -G "${OUT}/*.fastq.gz" >/dev/null; then
    local ARCH="${OUT}/.old_$(timestamp)"
    mkdir -p "$ARCH"
    mv "${OUT}"/*.fastq.gz "$ARCH"/ || true
    mv "${OUT}"/*manifest.json "$ARCH"/ 2>/dev/null || true
    mv "${OUT}"/chunk_checks.log "$ARCH"/ 2>/dev/null || true
    echo "[INFO] Archived existing FASTQs/manifests to: $ARCH"
  fi

  # checks log (TSV)
  CHK_LOG="${OUT}/chunk_checks.log"
  : > "$CHK_LOG"
  printf "#time\tgroup\tchunk\tattempt\tneed\tgot\ttotal\twant\n" >> "$CHK_LOG"

  # Build + filter class FASTAs (>=500 bp)
  concat_class_fasta prokaryote     Prokaryote     "$REFS" "$BASE"
  filter_fasta_minlen "${REFS}/Prokaryote.fasta" "${REFS}/Prokaryote.filtered.fasta" 500
  mv "${REFS}/Prokaryote.filtered.fasta" "${REFS}/Prokaryote.fasta"

  concat_class_fasta microeukaryote Microeukaryote "$REFS" "$BASE"
  filter_fasta_minlen "${REFS}/Microeukaryote.fasta" "${REFS}/Microeukaryote.filtered.fasta" 500
  mv "${REFS}/Microeukaryote.filtered.fasta" "${REFS}/Microeukaryote.fasta"

  concat_class_fasta virus          Virus          "$REFS" "$BASE"
  filter_fasta_minlen "${REFS}/Virus.fasta" "${REFS}/Virus.filtered.fasta" 500
  mv "${REFS}/Virus.filtered.fasta" "${REFS}/Virus.fasta"

  concat_class_fasta plasmid        Plasmid        "$REFS" "$BASE"
  filter_fasta_minlen "${REFS}/Plasmid.fasta" "${REFS}/Plasmid.filtered.fasta" 500
  mv "${REFS}/Plasmid.filtered.fasta" "${REFS}/Plasmid.fasta"

  # Per-record weights on filtered FASTAs
  make_perrec Prokaryote     "$REFS" "${ABUND}/prokaryote_generic.tsv"     "${PERREC}/prokaryote_generic_perrecord.tsv"
  make_perrec Microeukaryote "$REFS" "${ABUND}/microeukaryote_generic.tsv" "${PERREC}/microeukaryote_generic_perrecord.tsv"
  make_perrec Virus          "$REFS" "${ABUND}/virus_generic.tsv"          "${PERREC}/virus_generic_perrecord.tsv"
  make_perrec Plasmid        "$REFS" "${ABUND}/plasmid_generic.tsv"        "${PERREC}/plasmid_generic_perrecord.tsv"

  # Class pair counts (GENERIC)
  local N_PROK=$(( TOTAL_PAIRS_GENERIC * 56 / 100 ))
  local N_EUK=$((  TOTAL_PAIRS_GENERIC * 24 / 100 ))
  local N_VIR=$((  TOTAL_PAIRS_GENERIC * 10 / 100 ))
  local N_PLA=$((  TOTAL_PAIRS_GENERIC * 10 / 100 ))

  read -r PROK_R1 PROK_R2 < <(gen_group_chunked "Prokaryote"     "${REFS}/Prokaryote.fasta"     "${PERREC}/prokaryote_generic_perrecord.tsv"     "$N_PROK" "${TMP}/prok" "$TMP")
  read -r EUK_R1  EUK_R2  < <(gen_group_chunked "Microeukaryote" "${REFS}/Microeukaryote.fasta" "${PERREC}/microeukaryote_generic_perrecord.tsv" "$N_EUK" "${TMP}/euk" "$TMP")
  read -r VIR_R1  VIR_R2  < <(gen_group_chunked "Virus"          "${REFS}/Virus.fasta"          "${PERREC}/virus_generic_perrecord.tsv"          "$N_VIR" "${TMP}/vir" "$TMP")
  read -r PLA_R1  PLA_R2  < <(gen_group_chunked "Plasmid"        "${REFS}/Plasmid.fasta"        "${PERREC}/plasmid_generic_perrecord.tsv"        "$N_PLA" "${TMP}/pla" "$TMP")

  # GENERIC merge
  local GEN_TMP_R1="${TMP}/generic_R1.fastq.gz"
  local GEN_TMP_R2="${TMP}/generic_R2.fastq.gz"
  rm -f "$GEN_TMP_R1" "$GEN_TMP_R2"
  gz_repack_concat "$GEN_TMP_R1" "$PROK_R1" "$EUK_R1" "$VIR_R1" "$PLA_R1"
  gz_repack_concat "$GEN_TMP_R2" "$PROK_R2" "$EUK_R2" "$VIR_R2" "$PLA_R2"
  assert_pairs "$GEN_TMP_R1" "$GEN_TMP_R2" "$TOTAL_PAIRS_GENERIC" "GENERIC (TMP)"

  local GEN_R1="${OUT}/generic_R1.fastq.gz"
  local GEN_R2="${OUT}/generic_R2.fastq.gz"
  rm -f "$GEN_R1" "$GEN_R2"
  atom_mv "$GEN_TMP_R1" "$GEN_R1"
  atom_mv "$GEN_TMP_R2" "$GEN_R2"
  manifest_json "GENERIC" "$GEN_R1" "$GEN_R2" > "${OUT}/generic_manifest.json"
  echo "[GENERIC] ${COMM}  R1 pairs=$(npairs "$GEN_R1")  R2 pairs=$(npairs "$GEN_R2")"

  if [[ "$SHUFFLE" == "yes" ]]; then
    shuffle_pairs "$GEN_R1" "$GEN_R2" "${OUT}/generic_shuf_R1.fastq.gz" "${OUT}/generic_shuf_R2.fastq.gz"
  fi

  # FILTERED targets
  local F_PROK=$(( TOTAL_PAIRS_FILTERED * 14 / 100 ))
  local F_EUK=$((  TOTAL_PAIRS_FILTERED *  6 / 100 ))
  local F_VIR=$((  TOTAL_PAIRS_FILTERED * 40 / 100 ))
  local F_PLA=$((  TOTAL_PAIRS_FILTERED * 40 / 100 ))

  read -r FPROK_R1 FPROK_R2 < <(gen_group_chunked "Prokaryote"     "${REFS}/Prokaryote.fasta"     "${PERREC}/prokaryote_generic_perrecord.tsv"     "$F_PROK" "${TMP}/f_prok" "$TMP")
  read -r FEUK_R1  FEUK_R2  < <(gen_group_chunked "Microeukaryote" "${REFS}/Microeukaryote.fasta" "${PERREC}/microeukaryote_generic_perrecord.tsv" "$F_EUK"  "${TMP}/f_euk" "$TMP")
  read -r FVIR_R1  FVIR_R2  < <(gen_group_chunked "Virus"          "${REFS}/Virus.fasta"          "${PERREC}/virus_generic_perrecord.tsv"          "$F_VIR"  "${TMP}/f_vir" "$TMP")
  read -r FPLA_R1  FPLA_R2  < <(gen_group_chunked "Plasmid"        "${REFS}/Plasmid.fasta"        "${PERREC}/plasmid_generic_perrecord.tsv"        "$F_PLA"  "${TMP}/f_pla" "$TMP")

  local FIL_TMP_R1="${TMP}/filtered_R1.fastq.gz"
  local FIL_TMP_R2="${TMP}/filtered_R2.fastq.gz"
  rm -f "$FIL_TMP_R1" "$FIL_TMP_R2"
  gz_repack_concat "$FIL_TMP_R1" "$FPROK_R1" "$FEUK_R1" "$FVIR_R1" "$FPLA_R1"
  gz_repack_concat "$FIL_TMP_R2" "$FPROK_R2" "$FEUK_R2" "$FVIR_R2" "$FPLA_R2"
  assert_pairs "$FIL_TMP_R1" "$FIL_TMP_R2" "$TOTAL_PAIRS_FILTERED" "FILTERED (TMP)"

  local FIL_R1="${OUT}/filtered_R1.fastq.gz"
  local FIL_R2="${OUT}/filtered_R2.fastq.gz"
  rm -f "$FIL_R1" "$FIL_R2"
  atom_mv "$FIL_TMP_R1" "$FIL_R1"
  atom_mv "$FIL_TMP_R2" "$FIL_R2"
  manifest_json "FILTERED" "$FIL_R1" "$FIL_R2" > "${OUT}/filtered_manifest.json"
  echo "[FILTERED] ${COMM}  R1 pairs=$(npairs "$FIL_R1")  R2 pairs=$(npairs "$FIL_R2")"

  if [[ "$SHUFFLE" == "yes" ]]; then
    shuffle_pairs "$FIL_R1" "$FIL_R2" "${OUT}/filtered_shuf_R1.fastq.gz" "${OUT}/filtered_shuf_R2.fastq.gz"
  fi

  echo "[OK] ${COMM} short-read simulations complete."
  echo "====== [END] $COMM ======"
}

# -------------------- main --------------------
ARG="${1:-ALL}"

if [[ "$ARG" == "ALL" ]]; then
  for COMM_DIR in "${ROOT}"/train_communities/TrainC_short_*; do
    [[ -d "$COMM_DIR" ]] || continue
    COMM="$(basename "$COMM_DIR")"
    process_one "$COMM"
  done
else
  process_one "$ARG"
fi
