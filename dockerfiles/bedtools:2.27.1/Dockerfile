# TAG gcr.io/dnastack-dropbox-157402/private/profyle/bedtools:2.27.1
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	g++ \
	make \
	libz-dev \
	python-minimal \
	wget

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && \
	tar -zxvf bedtools-2.27.1.tar.gz --directory /opt/ && \
	rm bedtools-2.27.1.tar.gz
RUN cd /opt/bedtools2/ && \
	make && \
	make install
