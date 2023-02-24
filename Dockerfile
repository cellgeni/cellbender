# Start from nvidia-docker image with drivers pre-installed to use a GPU
FROM nvidia/cuda:11.3.1-cudnn8-runtime-ubuntu20.04

ENV DEBIAN_FRONTEND=noninteractive 

SHELL ["/bin/bash", "-c"]

# Install OS packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
    curl ca-certificates git bzip2 procps python3 python3-dev python3-pip python-is-python3 \
 && rm -rf /var/lib/apt/lists/*

# Install PyTorch
RUN pip install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu113

# Install Cellbender from source
RUN git clone --depth 1 --branch master https://github.com/broadinstitute/CellBender.git /opt/cellbender \
  && sed -i "s/version='0.2.0'/version='0.2.1'/g" /opt/cellbender/setup.py \
  && pip install -e /opt/cellbender --no-cache-dir \
  && rm -rf /tmp/*

# Install OS packages for QC script
RUN apt-get update && apt-get install -y liblzma-dev libbz2-dev zlib1g libpng-dev libxml2-dev \
    gfortran-7 libglpk-dev libhdf5-dev libcurl4-openssl-dev img2pdf wget libreadline8


# Install R for QC script
ARG R_VERSION=4.1.3

RUN wget https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-${R_VERSION}_1_amd64.deb && \
    apt-get update -qq && \
    apt-get install -f -y ./r-${R_VERSION}_1_amd64.deb && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/${R_VERSION}/lib/R /usr/lib/R && \
    rm r-${R_VERSION}_1_amd64.deb && \
    rm -rf /var/lib/apt/lists/*


# Install R packages for QC script
RUN Rscript -e "install.packages(c('argparse', 'stringr', 'hdf5r', 'Matrix'), repos='https://cloud.r-project.org/')"
