# TAG gcr.io/dnastack-dropbox-157402/private/profyle/annovar:2013.08.23
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	make \
	perl-modules-5.22 \
	locales

RUN locale-gen en_US.UTF-8
ENV LC_ALL en_US.UTF-8

COPY ./scripts/* /usr/local/bin/
