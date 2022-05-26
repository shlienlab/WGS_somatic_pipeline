# TAG gcr.io/dnastack-dropbox-157402/private/profyle/vcftools-perl:5.22.1
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

RUN apt-get -qq update && \
	apt-get -qq install \
		wget \
		autoconf \
		build-essential \
		pkg-config \
		zlib1g-dev \
		locales

RUN locale-gen en_US.UTF-8
ENV LC_ALL en_US.UTF-8

RUN wget https://github.com/vcftools/vcftools/archive/v0.1.16.tar.gz && \
	tar -zxvf v0.1.16.tar.gz --directory /opt/ && \
	rm v0.1.16.tar.gz

ENV PERL5LIB $PERL5LIB:/opt/vcftools-0.1.16/src/perl/

RUN cd /opt/vcftools-0.1.16 && \
	./autogen.sh && \
	./configure && \
	make && \
	make install

COPY ./scripts/* /usr/local/bin/
