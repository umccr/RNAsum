#!/bin/bash
set -euxo pipefail

echo "Run WTS-report"
Rscript /rmd_files/RNAseq_report.R --sample_name MDX190102_RNA010943 --dataset paad  --bcbio_rnaseq /work/WTS_data/MDX190102_RNA010943 --report_dir /work/output --umccrise /work/umccrise/PM3056445__MDX190101_DNA052297-T/ --ref_data_dir /work/WTS_ref_data/
