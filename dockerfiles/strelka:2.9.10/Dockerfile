# TAG gcr.io/dnastack-dropbox-157402/private/profyle/strelka:2.9.10
FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		wget \
		bzip2 \
		python

RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 && \
	tar -xjvf strelka-2.9.10.centos6_x86_64.tar.bz2 --directory /opt/ && \
	rm strelka-2.9.10.centos6_x86_64.tar.bz2
