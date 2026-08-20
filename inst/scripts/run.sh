#!/usr/bin/env bash

set -euo pipefail

# Example: run RNAsum in Docker (WGS + WTS) with the TEST reference dataset.
#
# Usage:
#   VERSION=2.0.4 OUTDIR=/path/to/output TESTDATA=/path/to/inst/rawdata/test_data ./run.sh # bundled test data
#   ---OR---
#   VERSION=2.0.4 OUTDIR=/path/to/output TESTDATA=/path/to/inputs ./run.sh                 # your own data
#
# Environment variables:
#   VERSION   RNAsum image tag to pull            (required)
#   OUTDIR    Host directory for the HTML report (required)
#   TESTDATA  Host directory with input files    (required)
#
# Input files are mounted read-only at /inputs; the report is written to
# ${OUTDIR}. Check https://github.com/umccr/RNAsum/tags for the latest tag.

if [[ -z "${VERSION:-}" ]]; then
  echo "Error: set VERSION to the RNAsum image tag to run." >&2
  echo "  See https://github.com/umccr/RNAsum/tags for available tags." >&2
  echo "  e.g. VERSION=2.0.4 OUTDIR=/path/to/output ./run.sh" >&2
  exit 1
fi
IMAGE="ghcr.io/umccr/rnasum:${VERSION}"

# Host directory holding the input data. Required: point it at your data, or at
# the bundled test data in this repo (inst/rawdata/test_data).
if [[ -z "${TESTDATA:-}" ]]; then
  echo "Error: set TESTDATA to the directory holding the input files." >&2
  echo "  e.g. TESTDATA=/path/to/inst/rawdata/test_data" >&2
  exit 1
fi

# Host directory for the generated HTML report. Required: set it explicitly so
# you know where the report lands.
if [[ -z "${OUTDIR:-}" ]]; then
  echo "Error: set OUTDIR to the directory where the report should be written." >&2
  echo "  e.g. OUTDIR=/path/to/output ./run.sh" >&2
  exit 1
fi
mkdir -p "${OUTDIR}"

docker pull "${IMAGE}"

# Mount input data at /inputs and outputs at /outputs.
# WGS + WTS using the TEST reference dataset.
docker run --rm \
  -v "${TESTDATA}:/inputs:ro" \
  -v "${OUTDIR}:/outputs" \
  "${IMAGE}" \
  rnasum.R \
    --sample_name test_sample_WTS \
    --dataset TEST \
    --salmon         "/inputs/dragen/TEST.quant.genes.sf" \
    --arriba_pdf     "/inputs/dragen/arriba/fusions.pdf" \
    --arriba_tsv     "/inputs/dragen/arriba/fusions.tsv" \
    --dragen_fusions "/inputs/dragen/test_sample_WTS.fusion_candidates.final" \
    --pcgr_tiers_tsv "/inputs/small_variants/TEST-snvs_indels.tiers.tsv" \
    --cn_gene_tsv    "/inputs/copy_number/TEST.cnv.gene.tsv" \
    --sv_tsv         "/inputs/structural/TEST-sv.tsv" \
    --report_dir /outputs

echo "Report written to ${OUTDIR}/"
