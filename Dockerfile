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
ARG LOCK_FILE="rnasum-linux-64.lock"
COPY ./deploy/conda/env/lock/${LOCK_FILE} .
RUN conda-lock install --name ${ENV_NAME} ${LOCK_FILE} && \
    mamba clean --all --force-pkgs-dirs

ARG MAMBA_PREFIX="/opt/conda"
ENV PATH="${MAMBA_PREFIX}/envs/${ENV_NAME}/bin:${PATH}"
ENV CONDA_PREFIX="${MAMBA_PREFIX}/envs/${ENV_NAME}"

CMD [ "rnasum.R" ]
