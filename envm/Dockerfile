FROM ubuntu:20.04
MAINTAINER Sehrish Kanwal "https://github.com/skanwal"

# Set up files
COPY ./envm/environment.yaml /tmp/environment.yaml
COPY rmd_files /rmd_files

# Runtime args
ARG TARGETPLATFORM
ARG BUILDPLATFORM

# Set args
ARG USER="rnasum_user"
ARG UID=1000
ARG GID=1000
ARG GROUP="rnasum_group"

# Set versions
ARG MAMBAFORGE_VERSION="4.10.3-4"
ARG PYPY_VERSION="3"
ARG PANDOC_VERSION="2.14.1"

ARG CRAN_KEY_SERVER="keyserver.ubuntu.com"
ARG CRAN_RECV_KEY="E298A3A825C0D65DFD57CBB651716619E084DAB9"

# Set env vars
ENV PATH /miniconda/bin:$PATH
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

# RUN through
RUN echo "Setting 'platform_arg'" 1>&2 && \
      echo "Got target platform value ${TARGETPLATFORM}" 1>&2 && \
      if [ "$TARGETPLATFORM" = "linux/amd64" ]; then \
        mamba_aws_platform_arg="x86_64"; \
        pandoc_platform_arg="amd64"; \
      elif [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
        mamba_aws_platform_arg="aarch64"; \
        pandoc_platform_arg="arm64"; \
      else \
        echo "Warning target platform $TARGETPLATFORM is not supported" 1>&2 && \
        echo "Defaulting to linux/amd64" 1>&2 && \
        mamba_aws_platform_arg="x86_64"; \
        pandoc_platform_arg="amd64"; \
      fi && \
      echo "Got platform arg value ${mamba_aws_platform_arg}" 1>&2

RUN echo "Updating and installing basics" 1>&2 && \
      DEBIAN_FRONTEND=noninteractive apt-get update -y -q && \
      DEBIAN_FRONTEND=noninteractive apt-get install -y -q \
        curl \
        wget \
        git \
        unzip \
        tar \
        gzip \
        bzip2 \
        g++ \
        make \
        zlib1g-dev \
        dbus \
        apt-transport-https \
        apt-utils \
        libcurl4-openssl-dev \
        libssl-dev \
        pkg-config \
        locales language-pack-en \
        tzdata && \
      echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

RUN echo "Install fonts for pandoc/rmarkdown - credit umccrise" 1>&2 && \
      echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections && \
        DEBIAN_FRONTEND=noninteractive apt-get update -y -q && \
        DEBIAN_FRONTEND=noninteractive apt-get install -y -q \
          ttf-mscorefonts-installer && \
    echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

RUN echo "Update locales" 1>&2 && \
      DEBIAN_FRONTEND=noninteractive apt-get update -y -q && \
      DEBIAN_FRONTEND=noninteractive locale-gen en_US.UTF-8 && \
      DEBIAN_FRONTEND=noninteractive dpkg-reconfigure locales && \
      echo "UTC" > /etc/timezone && \
      DEBIAN_FRONTEND=noninteractive dpkg-reconfigure tzdata

RUN echo "Install aws v2" 1>&2 && \
      curl "https://awscli.amazonaws.com/awscli-exe-linux-${mamba_aws_platform_arg}.zip" -o "awscliv2.zip" && \
      unzip awscliv2.zip && \
      ./aws/install && \
      rm -rf "aws/" "awscliv2.zip" \


RUN echo "Installing R dependencies" 1>&2 && \
      DEBIAN_FRONTEND=noninteractive apt-get install -y -q \
        dirmngr \
        gnupg \
        apt-transport-https \
        ca-certificates \
        software-properties-common && \
    echo "Installing Rpdftools dependencies" 1>&2 && \
        DEBIAN_FRONTEND=noninteractive apt-get install -y -q \
          libpoppler-cpp-dev && \
    echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

RUN echo "Add CRAN to sources list repository" 1>&2 && \
        DEBIAN_FRONTEND=noninteractive apt-key adv \
          --keyserver "${CRAN_KEY_SERVER}" \
          --recv-keys "${CRAN_RECV_KEY}" && \
        DEBIAN_FRONTEND=noninteractive add-apt-repository \
          'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

RUN echo "Installing pandoc" 1>&2 && \
        DEBIAN_FRONTEND=noninteractive apt-get update -y -q && \
        DEBIAN_FRONTEND=noninteractive apt-get install -y -q \
          texlive && \
        wget --no-verbose --output-document "pandoc.deb" \
        "https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-1-${pandoc_platform_arg}.deb" && \
        DEBIAN_FRONTEND=noninteractive dpkg --install pandoc.deb && \
        rm pandoc.deb && \
    echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

RUN echo "Install R" 1>&2 && \
        DEBIAN_FRONTEND=noninteractive apt-get install -y -q \
          r-base && \
    echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

RUN echo "Download and Install conda" 1>&2 && \
      wget --no-verbose --output-document mambaforge.sh \
        "https://github.com/conda-forge/miniforge/releases/download/${MAMBAFORGE_VERSION}/Mambaforge-pypy${PYPY_VERSION}-${MAMBAFORGE_VERSION}-Linux-${mamba_aws_platform_arg}.sh" && \
      bash mambaforge.sh -b -p /miniconda && \
    echo "Deleting conda install script" 1>&2 && \
      rm mambaforge.sh && \
    echo "Update conda" 1>&2 && \
      conda update --yes \
        --name base \
        --channel defaults \
        conda && \
    echo "Install environment with mamba" 1>&2 && \
      mamba env create -f /tmp/environment.yaml && \
    echo "Cleaning up conda" && \
      conda clean --all -y

RUN echo "Install R packages though R" 1>&2 && \
      Rscript -e 'options(Ncpus=2L, repos="https://cloud.r-project.org/"); \
                  install.packages(c( "BiocManager", "tidyverse", "rapportools", "optparse", \
                                      "openxlsx", "readr", "dplyr", "tidyr",\
                                      "glue", "rlang", "DT", "matrixStats", "knitr", \
                                      "plotly", "scales", "RCircos", \
                                      "ggplot2", "ggforce", "pdftools", \
                                      "png", "htmltools", "htmlwidgets", \
                                      "devtools", "manhattanly", "lares" ), \
                                   quiet=TRUE)' && \
    echo "Installing Bioconda packages through R" 1>&2 && \
      Rscript -e 'options(Ncpus=2L, repos="https://cloud.r-project.org/"); \
                  BiocManager::install(ask=FALSE, pkgs=c( "edgeR", "limma", "EDASeq", \
                                                          "preprocessCore", "tximport", "rhdf5", \
                                                          "EnsDb.Hsapiens.v75", \
                                                          "EnsDb.Hsapiens.v86", \
                                                          "BSgenome.Hsapiens.UCSC.hg19", \
                                                          "BSgenome.Hsapiens.UCSC.hg38" ), \
                                       quiet=TRUE)' && \
    echo "Clean up R downloaded packages" 1>&2 && \
    rm -rf "/tmp/downloaded_packages/" "/tmp/*.rds"

RUN echo "Adding user $USER:$GROUP" 1>&2 && \
      addgroup \
        --system \
        --gid "${GID}" \
        "${GROUP}" && \
      adduser \
        --system \
        --disabled-password \
        --ingroup "${GROUP}" \
        --uid "${UID}" \
        "$USER" && \
    echo "Changing ownership to user for /miniconda and /rmd_files" 1>&2 && \
      chown -R "$USER:$GROUP" /miniconda/ && \
      chown -R "$USER:$GROUP" /rmd_files/

RUN echo "Docker RUN time build complete!" 1>&2

# Set user
USER "$USER"

# Finalise paths for conda env
ENV PATH /miniconda/envs/rnasum/bin:/miniconda/bin:$PATH
ENV CONDA_PREFIX /miniconda/envs/rnasum

# Set runtime working directory
WORKDIR /rmd_files
