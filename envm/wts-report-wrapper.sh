#!/bin/bash
set -euxo pipefail

echo "Run WTS-report"
Rscript /rmd_files/RNAseq_report.R --sample_name test_sample_WTS  --bcbio_rnaseq /work/data/test_data/final/test_sample_WTS --report_dir /work/output --umccrise /work/data/test_data/umccrised/test_subject__test_sample_WGS --ref_data_dir /work/data --dataset TEST --save_tables FALSE
