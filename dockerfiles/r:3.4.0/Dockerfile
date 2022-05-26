# TAG gcr.io/dnastack-dropbox-157402/private/profyle/r:3.4.0
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

MAINTAINER Heather Ward <heather@dnastack.com>

# R 3.4.0
RUN echo "deb http://cloud.r-project.org/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN apt-get -qq update && \
	apt-get -qq install \
	r-base \
	r-base-dev

RUN R -e "install.packages('getopt')"
