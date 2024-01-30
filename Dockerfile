# syntax=docker/dockerfile:1

FROM ubuntu:jammy-20240111

# Copy code and libraries to docker image
RUN mkdir -p /home/lac
COPY code /home/lac

# Install compiler and required build tools
RUN apt-get update
RUN apt-get install -y gcc g++ m4 autoconf autotools-dev make libtool cmake

# Setup GMP library
RUN mkdir -p /home/lac/libs/gmp
ADD https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz /home/lac/libs/gmp/gmp.tar.xz
WORKDIR /home/lac/libs/gmp
RUN tar -xf gmp.tar.xz
RUN ./configure
RUN make
RUN make check
RUN make install

# Setup MPFR library
RUN mkdir -p /home/lac/libs/mpfr
ADD https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz /home/lac/libs/mpfr/mpfr.tar.xz
WORKDIR /home/lac/libs/mpfr
RUN tar -xf mpfr.tar.xz
RUN ./configure
RUN make
RUN make check
RUN make install

# Setup FLINT library
WORKDIR /home/lac/libs/flint
RUN git checkout flint-3.0
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
