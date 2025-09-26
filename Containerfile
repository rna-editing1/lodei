FROM docker.io/condaforge/mambaforge:latest 

RUN mkdir /lodei

COPY lodei /lodei/lodei
COPY pyproject.toml /lodei

#RUN mamba install -y python=3.10 && \
RUN mamba install -y -c bioconda pysamstats=1.1.2 samtools && \
    mamba install -y matplotlib pandas

RUN cd /lodei && \
    python -m pip install . --no-deps --ignore-installed

