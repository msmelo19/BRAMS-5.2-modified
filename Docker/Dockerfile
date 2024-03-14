# Build stage with Spack pre-installed and ready to be used
FROM spack/ubuntu-focal:latest as builder


# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - hdf5 +fortran +hl +mpi" \
&&   echo "  view: /opt/view" \
&&   echo "  concretizer:" \
&&   echo "    unify: false" \
&&   echo "  config:" \
&&   echo "    install_missing_compilers: true" \
&&   echo "    install_tree: /opt/software" ) > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && \
    spack env activate . && \
    spack install --fail-fast && \
    spack gc -y

# Strip all the binaries
RUN find -L /opt/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip -s

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

# Bare OS image to run the installed executables
FROM ubuntu:20.04

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/._view /opt/._view
COPY --from=builder /opt/view /opt/view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh

ARG DEBIAN_FRONTEND=noninteractive
ARG TZ=America/Sao_Paulo

RUN apt-get --yes -qq update \
    && apt-get --yes -qq install \
       build-essential \
       zlib1g-dev \
       gfortran \
       libopenmpi-dev \
       openmpi-bin \
       git \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/robertopsouto/BRAMS-5.2-modified.git
RUN cd BRAMS-5.2-modified && git checkout fixbugs

COPY comp.gnu.nojules.sh /BRAMS-5.2-modified/build_release
RUN cd /BRAMS-5.2-modified/build_release && ./comp.gnu.nojules.sh && make

COPY run-brams.sh /

ENTRYPOINT ["/bin/bash","./run-brams.sh"]

