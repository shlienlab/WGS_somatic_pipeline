# TAG gcr.io/dnastack-dropbox-157402/private/profyle/fastqc:0.11.5
FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	wget \
	unzip \
	build-essential \
	default-jre

# Install fastqc v 0.11.5
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
        unzip -d /opt fastqc_v0.11.5.zip && \
        rm fastqc_v0.11.5.zip && \
        chmod 777 /opt/FastQC/fastqc
RUN ln -s /opt/FastQC/fastqc /usr/local/bin/
