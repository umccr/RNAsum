FROM condaforge/mambaforge:23.1.0-4 as conda
LABEL maintainer="https://github.com/pdiakumis"

# install conda-lock
RUN mamba config \
      --set always_yes yes \
      --set always_copy yes && \
    mamba install \
      -c conda-forge \
      -c nodefaults \
      conda-lock && \
    mamba clean --all --force-pkgs-dirs

ARG ENV_NAME="rnasum_env"
COPY ./deploy/conda/env/lock/conda-lock.yml .
RUN conda-lock install --name ${ENV_NAME} conda-lock.yml && \
    mamba clean --all --force-pkgs-dirs && \
    rm conda-lock.yml

ARG MAMBA_PREFIX="/opt/conda"
ENV PATH="${MAMBA_PREFIX}/envs/${ENV_NAME}/bin:${PATH}"
ENV CONDA_PREFIX="${MAMBA_PREFIX}/envs/${ENV_NAME}"

CMD [ "rnasum.R" ]