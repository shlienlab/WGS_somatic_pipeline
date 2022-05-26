# TAG gcr.io/dnastack-dropbox-157402/private/profyle/delly:0.7.1
FROM ubuntu:14.04

MAINTAINER Heather Ward <heather@dnastack.com>

# install required packages
RUN apt-get -qq update && apt-get -qq install \
    build-essential \
    g++ \
    git \
    cmake \
    zlib1g-dev \
    ant \
    libbz2-dev \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# environment
ENV BOOST_ROOT /usr
ENV SEQTK_ROOT /htslib

# htslib (latest; last tested on 1.9)
RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib \
    && make \
    && make lib-static \
    && make install
ENV LD_LIBRARY_PATH /opt/htslib/

# delly 0.7.1
RUN wget https://github.com/dellytools/delly/archive/v0.7.1.tar.gz && \
    tar -zxvf v0.7.1.tar.gz --directory /opt/ && \
    rm v0.7.1.tar.gz
RUN cd /opt/delly-0.7.1 && \
    touch .htslib .boost && \
    make src/delly && \
    make src/stats

RUN ln -s /opt/delly-0.7.1/src/delly /usr/local/bin/