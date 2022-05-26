# TAG gcr.io/dnastack-dropbox-157402/private/profyle/gatk:3.8
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	openjdk-8-jre

# Version 3.8-0-ge9d806836
COPY ./GenomeAnalysisTK.jar /opt/gatk/3.8-0-ge9d806836/

ENV GATK /opt/gatk/3.8-0-ge9d806836/GenomeAnalysisTK.jar

