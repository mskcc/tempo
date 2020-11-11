# adapted from https://bitbucket.org/weischenfeldt/svaba/src/master/Dockerfile
FROM centos:7.3.1611

RUN yum install -y epel-release
RUN yum install -y \
    git \
    zlib-devel \
    gcc gcc-c++ \
    gcc-gfortran \
    make \
    wget \
    bzip2-devel.x86_64 \
    xz-devel

# they use 4a0606eba1bfd43c3d38433b27ac7f5e4887bc1e but that is outdated at this point so start with latest commit
# Submodule path 'SeqLib': checked out 'f7a89a127409a3f52fdf725fa74e5438c68e48fb'
# Submodule path 'SeqLib/bwa': checked out 'c02766e3c34ac3f4af9842e20a54b7f9f4b36d0b'
# Submodule path 'SeqLib/fermi-lite': checked out '5bc90f8d70e2b66184eccbd223a3be714c914365'
# Submodule path 'SeqLib/htslib': checked out '49fdfbda20acbd73303df3c7fef84f2d972c5f8d'
RUN git clone  https://github.com/walaj/svaba && \
cd svaba && \
git checkout 4d7b571356661c4ad34893ca2a47399212e0632b && \
git submodule update --init --recursive && \
./configure && \
make && \
make install

ENV PATH=/svaba/bin/:$PATH
