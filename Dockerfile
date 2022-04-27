# Start from nvidia-docker image with drivers pre-installed to use a GPU
FROM nvcr.io/nvidia/cuda:11.3.1-cudnn8-runtime-ubuntu20.04

SHELL ["/bin/bash", "-c"]

# Install OS packages
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    curl ca-certificates git bzip2 procps python3 python3-dev python3-pip python-is-python3 \
 && rm -rf /var/lib/apt/lists/*

# Install PyTorch
RUN pip install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu113

# Install Cellbender from source
RUN git clone --depth 1 --branch master https://github.com/broadinstitute/CellBender.git /opt/cellbender \
  && pip install -e /opt/cellbender --no-cache-dir \
  && rm -rf /tmp/*

