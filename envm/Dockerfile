FROM continuumio/miniconda3

ADD ./environment.yaml /tmp/environment.yaml
RUN conda env create -f /tmp/environment.yaml

#Pull the environment name out of the environment.yaml
RUN echo "source activate $(sed -n '5p' environment.yaml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(sed -n '5p' environment.yaml | cut -d' ' -f2)/bin:$PATH
