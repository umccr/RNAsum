#get miniconda installer
FROM continuumio/miniconda3

#copy environment file for conda and use it for creating conda environment
COPY ./envm/environment.yaml /tmp/environment.yaml
RUN conda env create -f /tmp/environment.yaml

#pull the environment name out of the environment.yaml, activate the conda envrionment and add to path
RUN echo "source activate $(sed -n '5p' /tmp/environment.yaml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/rnaseq-report/bin:$PATH

#get WTS report's RMD and R params script
COPY rmd_files /rmd_files
RUN chmod +x /rmd_files/*

#set working directory
WORKDIR /rmd_files

#define entry point
ENTRYPOINT ["Rscript", "/rmd_files/RNAseq_report.R"]
