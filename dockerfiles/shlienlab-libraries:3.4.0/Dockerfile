# TAG gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0
FROM gcr.io/dnastack-dropbox-157402/private/profyle/r:3.4.0

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		libssl-dev \
		libcurl4-gnutls-dev \
		libxml2-dev \
		cmake \
		git \
		pkg-config

# From https://github.com/r-lib/gert/issues/107
RUN git clone --depth=1 -b v1.0.0 https://github.com/libgit2/libgit2.git ~/libgit2_src \
    && cd ~/libgit2_src \
    && cmake . -DBUILD_CLAR=OFF -DCMAKE_BUILD_TYPE=Release -DEMBED_SSH_PATH=~/libssh2_src -DCMAKE_INSTALL_PREFIX=~/libgit2 \
    && cmake --build . --target install \
    && cp -r ~/libgit2/* /usr/bin \
    && cp -r ~/libgit2/* /usr/local

RUN apt-get -qq install libgit2-dev
#RUN Rscript -e 'install.packages("devtools")'

RUN Rscript -e 'install.packages("plyr")' \
	-e 'install.packages("dplyr")' \
	-e 'install.packages("vcfR")' \
	-e 'install.packages("stringr")' \
	-e 'install.packages("devtools")' \
	-e 'install.packages("readr")' \
	-e 'install.packages("tidyr")'

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' \
	-e 'biocLite("Biostrings")' \
	-e 'biocLite("biomaRt")' \
	-e 'biocLite("hiAnnotator")' \
	-e 'biocLite("Rsamtools")'

RUN Rscript -e 'install.packages("remotes")' \
	-e 'remotes::install_github("shlienlab/cosmic.cancer.gene.census")'

RUN mkdir /opt/libraries
COPY ./libraries/* /opt/libraries/

RUN Rscript -e 'install.packages("/opt/libraries/shlienlab.core.annotation.tar.gz", type="source")'
RUN Rscript -e 'install.packages("/opt/libraries/shlienlab.core.snv.tar.gz", type="source")'
RUN Rscript -e 'install.packages("/opt/libraries/shlienlab.core.filter.tar.gz", type="source")'
RUN Rscript -e 'install.packages("/opt/libraries/shlienlab.core.ssm.tar.gz", type="source")'
RUN Rscript -e 'install.packages("/opt/libraries/BICseq.tar.gz", type="source")'

COPY ./scripts/* /opt/scripts/
COPY ./ssm/* /opt/scripts/
