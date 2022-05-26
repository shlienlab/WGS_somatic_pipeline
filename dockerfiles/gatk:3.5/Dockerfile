# TAG gcr.io/dnastack-dropbox-157402/private/profyle/gatk:3.5
FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	software-properties-common \
	wget \
	bzip2

RUN add-apt-repository ppa:openjdk-r/ppa
RUN apt-get -qq update && \
	apt-get -qq install \
	openjdk-7-jdk

RUN update-alternatives --set java /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java

RUN wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.5-0-g36282e4" && \
	tar -jxvf "auth?package=GATK-archive&version=3.5-0-g36282e4" --directory /opt/ && \
	rm "auth?package=GATK-archive&version=3.5-0-g36282e4" && \
	rm -rf /opt/resources

ENV GATK /opt/GenomeAnalysisTK.jar