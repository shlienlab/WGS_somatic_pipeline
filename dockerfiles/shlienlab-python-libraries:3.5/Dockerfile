# TAG gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-python-libraries:3.5
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

RUN apt-get -qq update && \
	apt-get -qq install \
		wget \
		build-essential \
		zlib1g-dev \
		libbz2-dev \
		libz-dev \
		libssl-dev \
		libncurses5-dev \
		libbz2-dev \
		liblzma-dev \
		emboss \
		python3-pip

RUN pip3 install --upgrade pip
RUN pip3 install PyVCF numpy pysam

RUN mkdir -p /opt/{scripts,libraries}
COPY scripts/* /opt/scripts/
COPY libraries/* /opt/libraries/

ENV PYTHONPATH $PYTHONPATH:/opt/libraries
