# TAG gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-filtering:3.6
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
		emboss

RUN wget https://www.python.org/ftp/python/3.6.3/Python-3.6.3.tgz && \
    tar -xvf Python-3.6.3.tgz --directory /opt/ && \
    rm Python-3.6.3.tgz
RUN cd /opt/Python-3.6.3/ && \
    ./configure && \
    make && \
    make install

RUN pip3 install --upgrade pip
RUN pip3 install \
	pysam==0.15.2 \
	pandas==0.22.0 \
	pyyaml==5.1.1 \
	Biopython==1.73 \
	networkx==2.2 \
	matplotlib==2.2.4

# Install samtools 1.8
RUN wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 && \
        tar -jxvf samtools-1.8.tar.bz2 --directory /opt/ && \
        rm samtools-1.8.tar.bz2
RUN cd /opt/samtools-1.8 && \
        ./configure --prefix=/usr/local/ && \
        make && \
        make install

# Leave Shlien command line unchanged
RUN echo "alias python=python3.6" >> ~/.bashrc

RUN mkdir /opt/filterpipeline /opt/scripts
COPY ./filterpipeline /opt/filterpipeline/
COPY ./ssm /opt/scripts

ENV CLOUDSDK_PYTHON /usr/bin/python3.5
