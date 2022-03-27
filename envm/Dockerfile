FROM docker.io/ubuntu:20.04
MAINTAINER Sehrish Kanwal "https://github.com/skanwal"

# Runtime args
ARG TARGETPLATFORM
ARG BUILDPLATFORM

# Set args
ARG USER="rnasum_user"
ARG UID=1000
ARG GID=1000
ARG GROUP="rnasum_group"

# Set versions
ARG MAMBAFORGE_VERSION="4.12.0-0"
ARG PYPY_VERSION="3"
ARG PANDOC_VERSION="2.17.1.1"
ARG R_VERSION="4"
ARG R_PATCH_VERSION="4.1.3"
ARG R_SRC_URL="https://cran.rstudio.com/src/base/R-$R_VERSION/R-$R_PATCH_VERSION.tar.gz"
ARG R_INSTALL_PATH="/opt/R/$R_VERSION/$R_PATCH_VERSION"

ARG CRAN_KEY_SERVER="keyserver.ubuntu.com"
ARG CRAN_RECV_KEY="E298A3A825C0D65DFD57CBB651716619E084DAB9"

# Set env vars
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

RUN echo "Updating and installing basics" 1>&2 && \
      DEBIAN_FRONTEND=noninteractive apt-get update -y -qq && \
      DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
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
        libxml2-dev \
        pkg-config \
        locales language-pack-en \
        tzdata && \
      echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

RUN echo "Install fonts for pandoc/rmarkdown - credit umccrise" 1>&2 && \
      echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections && \
        DEBIAN_FRONTEND=noninteractive apt-get update -y -qq && \
        DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
          ttf-mscorefonts-installer && \
    echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

RUN echo "Update locales" 1>&2 && \
      DEBIAN_FRONTEND=noninteractive apt-get update -y -qq && \
      DEBIAN_FRONTEND=noninteractive locale-gen en_US.UTF-8 && \
      DEBIAN_FRONTEND=noninteractive dpkg-reconfigure locales && \
      echo "UTC" > /etc/timezone && \
      DEBIAN_FRONTEND=noninteractive dpkg-reconfigure tzdata

RUN echo "Install aws v2" 1>&2 && \
      echo "Got target platform value ${TARGETPLATFORM}" 1>&2; \
        if [ "$TARGETPLATFORM" = "linux/amd64" ]; then \
          aws_platform_arg="x86_64"; \
        elif [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
          aws_platform_arg="aarch64"; \
        else \
          echo "Warning target platform $TARGETPLATFORM is not supported" 1>&2 && \
          echo "Defaulting to linux/amd64" 1>&2 && \
          aws_platform_arg="x86_64"; \
        fi; \
      wget --no-verbose --output-document "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-$aws_platform_arg.zip"  && \
      unzip -qq awscliv2.zip && \
      ./aws/install && \
      rm -rf "aws/" "awscliv2.zip"

# From https://community.rstudio.com/t/installing-r-4-0-3-on-raspberry-pi-4-with-ubuntu-20-10/89369
RUN echo "Installing dependencies to build R from source" 1>&2 && \
      DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        openjdk-11-jdk \
        dirmngr \
        gnupg \
        apt-transport-https \
        ca-certificates \
        software-properties-common \
        gfortran \
        libreadline6-dev \
        libx11-dev \
        libxt-dev \
        libpng-dev \
        libjpeg-dev \
        libcairo2-dev \
        xvfb \
        libbz2-dev \
        libzstd-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        texinfo \
        texlive \
        texlive-fonts-extra \
        screen \
        wget \
        libpcre2-dev && \
    echo "Installing Rpdftools dependencies" 1>&2 && \
        DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
          libpoppler-cpp-dev && \
    echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

# Docs here: https://docs.rstudio.com/resources/install-r-source/
RUN echo "Installing R" 1>&2 && \
    wget --quiet \
      --output-document "R-${R_PATCH_VERSION}.tar.gz" \
      "${R_SRC_URL}" && \
    tar -xf "R-${R_PATCH_VERSION}.tar.gz" && \
    ( \
      cd "R-${R_PATCH_VERSION}" && \
      ./configure \
        --prefix "${R_INSTALL_PATH}" \
        --enable-memory-profiling \
        --enable-R-shlib \
        --with-blas \
        --with-lapack && \
      make && \
      make install ) && \
    echo "Verifying installation" 1>&2 && \
      "${R_INSTALL_PATH}/bin/R" --version && \
    echo "Link installation to /usr/local/bin" 1>&2 && \
    ln -s "${R_INSTALL_PATH}/bin/R" /usr/local/bin/R && \
    ln -s "${R_INSTALL_PATH}/bin/Rscript" /usr/local/bin/Rscript && \
    echo "Clean up source" 1>&2 && \
      rm -rf "R-${R_PATCH_VERSION}" "R-${R_PATCH_VERSION}.tar.gz"

# TODO - can probably remove from conda env if installed here
RUN echo "Installing pandoc" 1>&2 && \
      echo "Got target platform value ${TARGETPLATFORM}" 1>&2 && \
      if [ "$TARGETPLATFORM" = "linux/amd64" ]; then \
        pandoc_platform_arg="amd64"; \
      elif [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
        pandoc_platform_arg="arm64"; \
      else \
        echo "Warning target platform $TARGETPLATFORM is not supported" 1>&2 && \
        echo "Defaulting to linux/amd64" 1>&2 && \
        pandoc_platform_arg="amd64"; \
      fi && \
      DEBIAN_FRONTEND=noninteractive apt-get update -y -qq && \
      DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        texlive && \
      wget --no-verbose --output-document "pandoc.deb" \
      "https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-1-$pandoc_platform_arg.deb" && \
      DEBIAN_FRONTEND=noninteractive dpkg --install "pandoc.deb" && \
      rm "pandoc.deb" && \
    echo "Clean up apt" 1>&2 && \
      rm -rf "/var/lib/apt/lists"

# Set up conda env file
COPY ./envm/environment.yaml /tmp/environment.yaml

RUN echo "Download and Install conda" 1>&2 && \
      echo "Got target platform value ${TARGETPLATFORM}" 1>&2 && \
        if [ "$TARGETPLATFORM" = "linux/amd64" ]; then \
          mamba_platform_arg="x86_64"; \
        elif [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
          mamba_platform_arg="aarch64"; \
        else \
          echo "Warning target platform $TARGETPLATFORM is not supported" 1>&2 && \
          echo "Defaulting to linux/amd64" 1>&2 && \
          mamba_platform_arg="x86_64"; \
        fi; \
      wget --no-verbose --output-document mambaforge.sh \
        "https://github.com/conda-forge/miniforge/releases/download/${MAMBAFORGE_VERSION}/Mambaforge-pypy${PYPY_VERSION}-${MAMBAFORGE_VERSION}-Linux-$mamba_platform_arg.sh" && \
      bash mambaforge.sh -b -p /miniforge && \
    echo "Deleting conda install script" 1>&2 && \
      rm mambaforge.sh && \
    echo "Update mamba" 1>&2 && \
      /miniforge/bin/mamba update --quiet --yes \
        --name base \
        --channel conda-forge \
        mamba && \
    echo "Install environment with mamba" 1>&2 && \
      /miniforge/bin/mamba env create --quiet \
        --file /tmp/environment.yaml && \
    echo "Cleaning up conda" && \
      /miniforge/bin/conda clean --all --quiet --yes


# Add in RMD files
# Add in RMD files
COPY rmd_files /rmd_files

# Install BiocManager and renv
RUN echo "Install Renv and BiocManager packages though R" 1>&2 && \
      Rscript -e 'options(Ncpus=2L, repos="https://cloud.r-project.org/"); \
                  install.packages(c( "BiocManager", "renv", "yaml"), \
                                   quiet=TRUE)'

# Biocmanager pacakges that need to be installed prior to renv being activated
# See https://github.com/rstudio/renv/issues/934 for more details
# We've also added in rmd_files/renv_required_libraries.R
# To make sure that libraries that aren't auto downloaded by renv are added to the list
RUN echo "Install Bioconductor packages thorugh R" 1>&2 && \
      Rscript -e 'options(repos=BiocManager::repositories()); \
                  BiocManager::install(ask=FALSE, pkgs=c( \
                                                         "EDASeq", \
                                                         "limma", \
                                                         "preprocessCore", \
                                                         "edgeR", \
                                                         "rhdf5", \
                                                         "tximport", \
                                                         "EnsDb.Hsapiens.v75", \
                                                         "EnsDb.Hsapiens.v86", \
                                                         "BSgenome.Hsapiens.UCSC.hg19", \
                                                         "BSgenome.Hsapiens.UCSC.hg38" ), \
                                       quiet=TRUE)' && \
    echo "Clean up R downloaded packages" 1>&2 && \
    rm -rf "/tmp/downloaded_packages/" "/tmp/*.rds"


# Set the RENV_PATHS_CACHE environment variable to a shared path
# See for more info https://rstudio.github.io/renv/articles/renv.html#cache
ENV RENV_PATHS_CACHE="${R_INSTALL_PATH}/renv/shared"
RUN mkdir -p "${RENV_PATHS_CACHE}" && \
    chmod -R 777 "${RENV_PATHS_CACHE}"

# Install through renv
# Trying restart
RUN ( \
      cd "/rmd_files" && \
      if [ -f "renv.lock" ]; then \
        echo "Restoring the environment from the lock file" 1>&2; \
        Rscript -e "options(repos=BiocManager::repositories()); \
                    renv::restore()"; \
      else \
        echo "No lock file exists, initialising the environment" 1>&2; \
        Rscript -e 'options(repos=BiocManager::repositories()); \
                    renv::init()'; \
      fi; \
    )  

# Add non-root user
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
    echo "Changing ownership to user for /miniforge and /rmd_files" 1>&2 && \
      chown -R "$USER:$GROUP" /miniforge/ && \
      chown -R "$USER:$GROUP" /rmd_files/

RUN echo "Docker RUN time build complete!" 1>&2

# Finalise paths for conda env (but we don't need it just yet)
# Set user
USER "$USER"

# Set env for user
ENV PATH="/miniforge/envs/rnasum/bin:/miniforge/bin:$PATH"
ENV CONDA_PREFIX="/miniforge/envs/rnasum"
ENV USER="$USER"

# Solves https://community.rstudio.com/t/timedatectl-had-status-1/72060/4
ENV TZ="UTC"

# Set /rmd_files as a volume
VOLUME /rmd_files

# Set /rmd_files as the default working directory
WORKDIR /rmd_files
