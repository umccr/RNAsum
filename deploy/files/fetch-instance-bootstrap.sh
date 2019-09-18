#!/bin/bash
set -euxo pipefail

STACK="rna-seq_report"
SCRIPT="bootstrap-instance.sh"

# Fetch the bootstrap script 
wget https://raw.githubusercontent.com/umccr/workflows/master/$STACK/$SCRIPT -O "/$SCRIPT" && chmod +x "/$SCRIPT" && /bin/bash "/$SCRIPT"
