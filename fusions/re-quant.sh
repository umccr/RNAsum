#!/bin/bash
#script for requantifying pizzly results
#recommend running this in an interactive session
#expects kallisto in the path (recommend conda installation "conda install -c bioconda kallisto")

#check input parameters
while getopts "t:f:1:2:" option
do
case "${option}"
in
t) TRANSCRIPTOME=${OPTARG};;
f) FUSIONSFASTA=${OPTARG};;
1) FASTQ1=${OPTARG};;
2) FASTQ2=$OPTARG;;
esac
done

if [[ $TRANSCRIPTOME == "" ]] || [[ $FUSIONSFASTA == "" ]] || [[ $FASTQ1 == "" ]] || [[ $FASTQ2 == "" ]]; then
    echo -e "Positional parameters missing\nUsage: Re_quant.sh -t /path_to/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.fa -f /path_to_bcbio_WTS_run/final/*/pizzly/*.fusions.fasta -1 /path_to/*)R1_001.fastq.gz -2 /path_to/*_R2_001.fastq.gz" 
fi

#create a new index based on the transcriptome and the fusion transctripts identified by pizzly
cat $TRANSCRIPTOME $FUSIONSFASTA | gzip - > transcripts_with_fusions_fasta.gz
kallisto index -k 31 -i ./transcripts_with_fusions.kidx transcripts_with_fusions_fasta.gz

#run kallisto in normal quantification mode on the expanded index to quantify both normal transcripts and fusions
kallisto quant -i ./transcripts_with_fusions.kidx -o ./quant_pizzly_post $FASTQ1 $FASTQ2 
