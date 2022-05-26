# TAG gcr.io/dnastack-dropbox-157402/manta:1.6.0
FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	wget \
	bzip2 \
	python

RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 && \
	tar -xjvf manta-1.6.0.centos6_x86_64.tar.bz2 --directory /opt/ && \
	rm manta-1.6.0.centos6_x86_64.tar.bz2
