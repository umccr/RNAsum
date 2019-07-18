#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

hints:
 DockerRequirement:
  dockerPull: quay.io/biocontainers/kallisto:0.44.0--h7d86c95_2

inputs:
  transcriptome: File
  fusionsfasta: File
  reads: File[]

steps:
  concat:
    run: ./tools/cat.cwl
    in:
      transcriptome: transcriptome
      fusionsfasta: fusionsfasta
    out:
      - reference

  zip:
    run: ./tools/zip.cwl
    in:
      transcriptsWithFusions: concat/reference
    out:
      - zippedReference

  indexing:
    run: ./tools/kallisto-index.cwl
    in:
      fasta: zip/zippedReference
    out: 
      - index
      
  quantifying:
    run: ./tools/kallisto-quant.cwl
    in:
      fastqs: reads
      index: indexing/index
    out:
      - quantification

outputs:
  quant:
    type: File
    outputSource: quantifying/quantification


      
