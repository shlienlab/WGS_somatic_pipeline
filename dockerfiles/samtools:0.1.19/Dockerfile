# TAG gcr.io/dnastack-dropbox-157402/private/profyle/samtools:0.1.19
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

RUN wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2 && \
	tar -jxvf samtools-0.1.19.tar.bz2 --directory /opt/ && \
	rm samtools-0.1.19.tar.bz2
RUN cd /opt/samtools-0.1.19 && \
	make CFLAGS='-g -Wall -O2 -fPIC'
RUN ln -s /opt/samtools-0.1.19/samtools /usr/local/bin/
