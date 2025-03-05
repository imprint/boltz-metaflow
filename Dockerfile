# The base image may need to change depending on what version of CUDA the host has.
# FROM nvidia/cuda:12.4.0-runtime-ubuntu22.04
FROM nvidia/cuda:12.2.2-cudnn8-runtime-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

# Install dependencies as root
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    wget \
    jq \
    libgomp1 \
    libxml2-dev \
    zip \
    default-jre \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a new user named 'user'
RUN useradd -m -s /bin/bash user

# Switch to the non-root user
USER user
WORKDIR /home/user

# Conda setup
RUN wget -nv https://repo.anaconda.com/miniconda/Miniconda3-py311_24.11.1-0-Linux-x86_64.sh \
    && bash Miniconda3-*.sh -b -p /home/user/miniconda3 \
    && rm -f Miniconda3-*.sh \
    && /home/user/miniconda3/bin/conda init bash \
    && /home/user/miniconda3/bin/conda config --add channels defaults \
    && /home/user/miniconda3/bin/conda config --add channels bioconda \
    && /home/user/miniconda3/bin/conda config --add channels conda-forge \
    && /home/user/miniconda3/bin/conda config --set channel_priority strict \
    && /home/user/miniconda3/bin/conda --version
ENV PATH="/home/user/miniconda3/bin:$PATH"

# Install boltz
RUN git clone --branch main https://github.com/timodonnell/boltz.git /home/user/boltz \
    && cd /home/user/boltz \
    && git remote add upstream https://github.com/jwohlwend/boltz.git \
    && pip install -e /home/user/boltz

# Test boltz (GPU unavailable error is expected but okay).
# We do this so it is forced to download the model weights.
RUN mkdir /tmp/test-data \
    && echo ">A|protein|empty" >> /tmp/test-data/protein.fasta \
    && echo "IVMTQSPATLSLSPGERATLSCRASQSAGFYLAWYQQKPGQAPRLLIYDTSNRATGIPARFSGRGSGTDFT" >> /tmp/test-data/protein.fasta \
    && boltz predict /tmp/test-data/protein.fasta ; ls -l /tmp/test-data

# Install additional python packages
RUN which pip \
    && pip install --upgrade pip setuptools wheel \
    && pip install --no-cache-dir \
      DockQ \
      ProDy

# Copy analysis files
COPY --chown=user . /home/user/analysis
WORKDIR /home/user/analysis
RUN mkdir /home/user/analysis/work
