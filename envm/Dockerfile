#get miniconda installer
FROM continuumio/miniconda3
MAINTAINER Sehrish Kanwal "https://github.com/skanwal"

#copy environment file for conda and use it for creating conda environment
COPY ./envm/environment.yaml /tmp/environment.yaml
RUN conda env create -f /tmp/environment.yaml

#pull the environment name out of the environment.yaml, activate the conda envrionment and add to path
RUN echo "source activate $(sed -n '5p' /tmp/environment.yaml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/rnaseq-report/bin:$PATH

#get WTS report's RMD and R params script
COPY rmd_files /rmd_files

#set working directory
WORKDIR /rmd_files

#copy wrapper scripts
COPY envm/wts-report-wrapper.sh /rmd_files/wts-report-wrapper.sh
RUN chmod +x /rmd_files/*

#set-up cmd for easier access
CMD ["/bin/bash", "wts-report-wrapper.sh"]
