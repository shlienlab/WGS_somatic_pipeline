# TAG gcr.io/dnastack-dropbox-157402/private/profyle/bwa-sambamba:0.7.8
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		wget \
		make \
		gcc \
		build-essential \
		libz-dev \
		libncurses5-dev \
		libbz2-dev \
		liblzma-dev

# BWA v0.7.8
RUN wget https://github.com/lh3/bwa/archive/0.7.8.tar.gz && \
	tar -zxvf 0.7.8.tar.gz --directory /opt/ && \
	rm 0.7.8.tar.gz
RUN cd /opt/bwa-0.7.8 && \
	make
RUN ln -s /opt/bwa-0.7.8/bwa /usr/local/bin/

# Sambamba v0.7.0
RUN wget -P /opt/ https://github.com/biod/sambamba/releases/download/v0.7.0/sambamba-0.7.0-linux-static.gz && \
	gunzip /opt/sambamba-0.7.0-linux-static.gz && \
	chmod +x /opt/sambamba-0.7.0-linux-static && \
	ln -s /opt/sambamba-0.7.0-linux-static /usr/local/bin/sambamba
