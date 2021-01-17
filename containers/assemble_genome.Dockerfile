FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - assemble_genome"
LABEL software.version="1.5.6"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/assemble_genome.yml"
LABEL conda.md5="9d01230a850123fc550d2ed0f600f04d"

COPY conda/linux/assemble_genome.yml /
RUN conda env create -q -f assemble_genome.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-assemble_genome/bin:$PATH