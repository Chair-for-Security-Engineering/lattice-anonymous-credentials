# syntax=docker/dockerfile:1

FROM ubuntu:jammy-20240111

# Copy code and libraries to docker image
RUN mkdir -p /home/lac
COPY code /home/lac

# Install compiler and required build tools
RUN apt-get update
RUN apt-get install -y gcc g++ m4 autoconf autotools-dev make libtool cmake xz-utils

# download and unpack GMP and MPFR
WORKDIR /home/lac/libs
ADD https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz /home/lac/libs/gmp.tar.xz
ADD https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz /home/lac/libs/mpfr.tar.xz
RUN tar -xf gmp.tar.xz
RUN tar -xf mpfr.tar.xz

# Setup GMP library
WORKDIR /home/lac/libs/gmp-6.2.1
RUN ./configure
RUN make
RUN make check
RUN make install

# Setup MPFR library
WORKDIR /home/lac/libs/mpfr-4.2.1
RUN ./configure
RUN make
RUN make check
RUN make install

# Setup FLINT library
WORKDIR /home/lac/libs/flint
RUN ./bootstrap.sh
RUN ./configure
RUN make -j
RUN make check
RUN make install

# Build lac
RUN mkdir -p /home/lac/_build
WORKDIR /home/lac/_build
RUN cmake ..
RUN make
