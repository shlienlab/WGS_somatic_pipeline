# TAG gcr.io/dnastack-dropbox-157402/private/profyle/bwa-samtools:0.7.8
FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		wget \
		make \
		gcc \
		build-essential \
		libz-dev \
		libncurses5-dev \
		libbz2-dev \
		liblzma-dev

# Bwa v0.7.8
RUN wget https://github.com/lh3/bwa/archive/0.7.8.tar.gz && \
	tar -zxvf 0.7.8.tar.gz --directory /opt/ && \
	rm 0.7.8.tar.gz
RUN cd /opt/bwa-0.7.8 && \
	make
RUN ln -s /opt/bwa-0.7.8/bwa /usr/local/bin/

# Samtools 1.5
RUN wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 && \
	tar -jxvf samtools-1.5.tar.bz2 --directory /opt/ && \
	rm samtools-1.5.tar.bz2
RUN cd /opt/samtools-1.5 && \
	./configure --prefix=/usr/local/ && \
	make && \
	make install