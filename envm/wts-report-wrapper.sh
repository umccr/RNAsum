#!/bin/bash
set -euxo pipefail

echo "Run WTS-report"
Rscript /rmd_files/RNAseq_report.R --sample_name test_sample_WTS  --dataset paad  --count_file ./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv --report_dir ./data/test_data/final/test_sample_WTS/RNAseq_report  --umccrise ./data/test_data/umccrised/test_sample_WGS --ref_data_dir /rmd_files/data