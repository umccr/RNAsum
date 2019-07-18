#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ "kallisto", "index" ]

hints:
 DockerRequirement:
  dockerPull: quay.io/biocontainers/kallisto:0.44.0--h7d86c95_2

requirements:
  ResourceRequirement:
    ramMin: 4096

inputs:
 fasta:
   type: File
   inputBinding: {}

arguments: [ --index, transcripts_with_fusions.kidx]

outputs:
 index:
   type: File
   outputBinding: 
     glob: transcripts_with_fusions.kidx
