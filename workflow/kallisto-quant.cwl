#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [ "kallisto", "quant" ]

hints:
 DockerRequirement:
  dockerPull: quay.io/biocontainers/kallisto:0.44.0--h7d86c95_2

requirements:
  ResourceRequirement:
    ramMin: 4096

inputs:
 fastqs:
   type: File[]
   inputBinding: {}

 index:
   type: File
   inputBinding:
     prefix: "--index"

arguments: [ "--output-dir", out ]
  
outputs:
 quantification:
  type: File
  outputBinding:
   glob: out/abundance.tsv
