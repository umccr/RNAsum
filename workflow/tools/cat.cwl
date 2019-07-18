#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: cat
stdout: transcripts_with_fusions.fasta
inputs:
  transcriptome:
    type: File
    inputBinding:
      position: 1
  fusionsFasta:
    type: File
    inputBinding:
      position: 2
outputs:
  reference:
    type: stdout
