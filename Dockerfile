FROM ubuntu:22.04

# Instala pacotes basicos
RUN set -ex && \
    apt-get update && \
    apt-get install -y lsb-release wget software-properties-common zlib1g-dev libomp-dev
  
# Instala cmake, ninja, gfortran e openblas
RUN set -ex && \
    apt-get update && \
    apt-get install -y cmake ninja-build gfortran libopenblas-dev