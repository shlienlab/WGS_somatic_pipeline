# TAG gcr.io/dnastack-dropbox-157402/private/profyle/samtools:1.9
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

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

RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2 && \
	tar -jxvf samtools-1.9.tar.bz2 --directory /opt/ && \
	rm samtools-1.9.tar.bz2
RUN cd /opt/samtools-1.9 && \
	make
	# make CFLAGS='-g -Wall -O2 -fPIC'
RUN ln -s /opt/samtools-1.9/samtools /usr/local/bin/
