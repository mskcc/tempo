# adapted from https://bitbucket.org/weischenfeldt/svaba/src/master/Dockerfile
# they use SvABA version 134 (4a0606eba1bfd43c3d38433b27ac7f5e4887bc1e) but that is outdated at this point so start with latest commit
FROM centos:7.3.1611

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.1"

RUN yum install -y epel-release
RUN yum install -y \
    git \
    zlib-devel \
    gcc gcc-c++ \
    gcc-gfortran \
    make \
    wget \
    bzip2-devel.x86_64 \
    xz-devel \
    bcftools

# original commit for pcawg
ENV SVABA_COMMIT 4a0606eba1bfd43c3d38433b27ac7f5e4887bc1e
# latest release commit
#ENV SVABA_COMMIT d12cf224f7a488b913eabbcf54a215e17238032c

# SvABA Version: 1.1.3 (4d7b57);
# Submodule path 'SeqLib': checked out 'f7a89a127409a3f52fdf725fa74e5438c68e48fb'
# Submodule path 'SeqLib/bwa': checked out 'c02766e3c34ac3f4af9842e20a54b7f9f4b36d0b'
# Submodule path 'SeqLib/fermi-lite': checked out '5bc90f8d70e2b66184eccbd223a3be714c914365'
# Submodule path 'SeqLib/htslib': checked out '49fdfbda20acbd73303df3c7fef84f2d972c5f8d'
RUN git clone  https://github.com/walaj/svaba && \
cd svaba && \
git checkout ${SVABA_COMMIT} && \
git submodule update --init --recursive && \
./configure && \
make && \
make install

ENV PATH=/svaba/bin/:$PATH
