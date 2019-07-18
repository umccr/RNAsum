#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [gzip, -c]
stdout: transcripts_with_fusions_fasta.gz
inputs:
  transcriptsWithFusions:
    type: File
    inputBinding:
      position: 1
outputs:
  zippedReference:
    type: stdout
