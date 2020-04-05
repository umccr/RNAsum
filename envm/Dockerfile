FROM ubuntu:16.04
MAINTAINER Sehrish Kanwal "https://github.com/skanwal"

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget git unzip tar gzip bzip2 g++ make zlib1g-dev dbus

# Install fonts for pandoc/rmarkdown - credit umccrise
RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections
RUN apt-get install -y ttf-mscorefonts-installer

# Install conda
RUN wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /miniconda
ENV PATH /miniconda/bin:$PATH

# Copy environment file for conda and use it for creating conda environment
COPY ./envm/environment.yaml /tmp/environment.yaml
RUN conda env create -f /tmp/environment.yaml

# Activate the environment and set path
ENV PATH  /miniconda/envs/rnasum/bin:$PATH
ENV CONDA_PREFIX /miniconda/envs/rnasum

#get WTS report's RMD and R params script
COPY rmd_files /rmd_files

#set working directory
WORKDIR /rmd_files

#change permissions
#RUN chmod +x /rmd_files/*
RUN chmod -R 777 /rmd_files

# Setting locales and timezones, based on https://github.com/jacksoncage/node-docker/blob/master/Dockerfile
# (setting UTC for readr expecting UTC https://rdrr.io/github/tidyverse/readr/src/R/locale.R - credit umccrise)
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
RUN apt-get update && \
    apt-get install -y locales language-pack-en && \
    locale-gen en_US.UTF-8 && \
    dpkg-reconfigure locales && \
    apt-get install -y tzdata && \
    echo "Etc/UTC" > /etc/timezone && \
    dpkg-reconfigure -f noninteractive tzdata
