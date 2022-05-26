# TAG gcr.io/dnastack-dropbox-157402/private/profyle/picard:1.108
FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	software-properties-common \
	wget \
	ant \
	subversion \
	unzip

RUN add-apt-repository ppa:openjdk-r/ppa
RUN apt-get -qq update && \
	apt-get -qq install \
	openjdk-7-jdk

RUN update-alternatives --set java /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java
ENV JAVA6_HOME /usr/lib/jvm/java-7-openjdk-amd64/

RUN mkdir -p /root/.ant/lib

# bcel required by ant to build Picard
RUN wget http://www.java2s.com/Code/JarDownload/bcel/bcel-5.2.jar.zip && \
	unzip bcel-5.2.jar.zip -d /root/.ant/lib/ && \
	rm bcel-5.2.jar.zip

# Picard tools 1.108
RUN wget https://github.com/broadinstitute/picard/archive/1.108.tar.gz && \
	tar -zxvf 1.108.tar.gz --directory /opt/ && \
	rm 1.108.tar.gz
RUN cd /opt/picard-1.108 && \
	ant

