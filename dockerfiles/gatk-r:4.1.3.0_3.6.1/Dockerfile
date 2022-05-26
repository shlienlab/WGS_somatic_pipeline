# TAG gcr.io/dnastack-dropbox-157402/private/profyle/gatk-r:4.1.3.0_3.6.1
FROM gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.3.0

RUN apt-get -qq update && \
	apt-get -qq install \
		make \
		gcc \
		build-essential \
		gfortran \
		libz-dev \
		libbz2-dev \
		liblzma-dev \
		libpcre3-dev \
		libcurl4-openssl-dev

# R 3.6.1
RUN wget https://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz && \
	tar -zxvf R-3.6.1.tar.gz --directory /opt/ && \
	rm R-3.6.1.tar.gz
RUN cd /opt/R-3.6.1 && \
	./configure --with-readline=no --with-x=no && \
	make && \
	make install
