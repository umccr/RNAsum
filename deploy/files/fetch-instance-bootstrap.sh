#!/bin/bash
set -euxo pipefail

# Install and start docker service
sudo yum update -y
sudo amazon-linux-extras install docker -y
sudo service docker start

# Add the ssm-user to the docker group 
sudo usermod -a -G docker ssm-user

STACK="rna-seq_report"
SCRIPT="bootstrap-instance.sh"

# Fetch the bootstrap script 
wget https://raw.githubusercontent.com/umccr/workflows/master/$STACK/$SCRIPT -O "/$SCRIPT" && chmod +x "/$SCRIPT"