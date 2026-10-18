#!/usr/bin/env bash
set -euo pipefail

projpath="/media/bafj7/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig"
indir="${projpath}"
outdir="${projpath}/derived_bigwigs"

p=0.001
binsize=50
max_jobs=8

mkdir -p "$outdir"

run_job () {
  local treat="$1" cond="$2" rep="$3"
  local L="${indir}/LB1_${treat}_${cond}_${rep}.bigWig"
  local M="${indir}/MS_${treat}_${cond}_${rep}.bigWig"
  local out="${outdir}/LB1_${treat}_${cond}_${rep}.adj_log2.bw"

  if [[ -f "$L" && -f "$M" ]]; then
    echo "[$(date +%H:%M:%S)] $treat $cond $rep"
    bigwigCompare -b1 "$L" -b2 "$M" \
      --operation log2 --pseudocount "$p" --binSize "$binsize" \
      -o "$out"
  else
    echo "Missing: $L or $M" >&2
  fi
}

for treat in DMSO ATMi; do
  for rep in R1 R2 R3; do
    for cond in CTRL LAD; do
      run_job "$treat" "$cond" "$rep" &

      # limiter: don’t exceed max_jobs running processes
      while :; do
        n_running=$(jobs -r | wc -l)
        if (( n_running < max_jobs )); then
          break
        fi
        sleep 0.2
      done
    done
  done
done

wait
echo "🎉 Done."
