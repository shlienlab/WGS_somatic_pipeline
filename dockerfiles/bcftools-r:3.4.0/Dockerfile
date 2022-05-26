# TAG gcr.io/dnastack-dropbox-157402/private/profyle/bcftools-r:3.4.0
FROM gcr.io/dnastack-dropbox-157402/private/profyle/r:3.4.0

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	wget

# Bcftools 1.6
RUN wget https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2 && \
	tar -jxvf bcftools-1.6.tar.bz2 --directory /opt/ && \
	rm bcftools-1.6.tar.bz2
RUN cd /opt/bcftools-1.6 && \
	./configure && \
	make && \
	make install