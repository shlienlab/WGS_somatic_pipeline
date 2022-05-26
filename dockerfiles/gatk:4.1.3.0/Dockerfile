# TAG gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.3.0
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	software-properties-common \
	wget \
	bzip2 \
	unzip \
	python \
	openjdk-8-jdk \
	tabix \
	build-essential \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl4-openssl-dev

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.3.0/gatk-4.1.3.0.zip && \
	unzip gatk-4.1.3.0.zip -d /opt/ && \
	rm gatk-4.1.3.0.zip
RUN ln -s /opt/gatk-4.1.3.0/gatk /usr/local/bin/

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
	tar -jxvf bcftools-1.9.tar.bz2 --directory /opt/ &&\
	rm bcftools-1.9.tar.bz2

RUN cd /opt/bcftools-1.9 && \
	./configure && \
	make && \
	make install
