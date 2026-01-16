FROM fedora:latest

ARG BUILD_TYPE=Release

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc-c++ \
        gdb \
        git \
        make \
        tree \
    && dnf clean all

# Copy the MIAM code
COPY . miam

# Build the library
RUN cd miam \
    && cmake -S . \
             -B build \
             -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
    && cd build \
    && make -j 5

WORKDIR miam/build