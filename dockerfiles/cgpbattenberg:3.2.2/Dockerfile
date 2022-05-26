# TAG gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2
FROM gcr.io/dnastack-dropbox-157402/private/profyle/r:3.4.0

RUN apt-get -qq update && \
	apt-get -qq install \
    build-essential \
    libbz2-dev \
    libcurl4-gnutls-dev \
    libgnutls-dev \
    liblzma-dev \
    libncurses5-dev \
    libssl-dev \
    locales \
    nettle-dev \
    zlib1g-dev \
    curl \
    libexpat1-dev \
    libgoogle-perftools-dev \
    lsof \
    time \
    psmisc \
    wget

# cgpBigWig 1.1.0
RUN wget https://github.com/cancerit/cgpBigWig/archive/1.1.0.tar.gz && \
	tar -zxvf 1.1.0.tar.gz --directory /opt/ && \
	rm 1.1.0.tar.gz
RUN cd /opt/cgpBigWig-1.1.0 && \
	./setup.sh /usr/local/

# PCAP-core 4.3.5
RUN wget https://github.com/cancerit/PCAP-core/archive/4.3.5.tar.gz && \
	tar -zxvf 4.3.5.tar.gz --directory /opt/ && \
	rm 4.3.5.tar.gz
RUN cd /opt/PCAP-core-4.3.5 && \
	./setup.sh /usr/local/

ENV PATH /usr/local/biobambam2/bin:$PATH
ENV PERL5LIB /usr/local/lib/perl5:$PERL5LIB

# alleleCount 4.0.0
RUN wget https://github.com/cancerit/alleleCount/archive/v4.0.0.tar.gz && \
	tar -zxvf v4.0.0.tar.gz --directory /opt/ && \
	rm v4.0.0.tar.gz
RUN cd /opt/alleleCount-4.0.0 && \
	./setup.sh /usr/local/

# cgpVcf 2.0.4
RUN wget https://github.com/cancerit/cgpVcf/archive/v2.0.4.tar.gz && \
	tar -zxvf v2.0.4.tar.gz --directory /opt/ && \
	rm v2.0.4.tar.gz
RUN cd /opt/cgpVcf-2.0.4 && \
	./setup.sh /usr/local/

# cgpBattenberg 3.2.2
RUN wget https://github.com/cancerit/cgpBattenberg/archive/v3.2.2.tar.gz && \
	tar -zxvf v3.2.2.tar.gz --directory /opt/ && \
	rm v3.2.2.tar.gz
RUN cd /opt/cgpBattenberg-3.2.2 && \
	./setup.sh /usr/local/
RUN cd /opt/cgpBattenberg-3.2.2/Rsupport && \
	./setupR.sh /usr/local/lib/R/site-library

# Impute2 2.3.1
RUN wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.1_x86_64_static.tgz && \
	tar -zxvf impute_v2.3.1_x86_64_static.tgz --directory /opt/ && \
	rm impute_v2.3.1_x86_64_static.tgz
RUN ln -s /opt/impute_v2.3.1_x86_64_static/impute2 /usr/local/bin/

# .Rprofile
RUN echo "options(bitmapType = 'cairo')" >> /etc/R/.Rprofile && \
	echo "options(scipen = 999)" >> /etc/R/.Rprofile

# .Renviron - cgpBattenberg installs in an odd location
RUN echo "R_LIBS_USER=/usr/local/lib/R/site-library/R-lib" >> /etc/R/.Renviron

ENV R_LIBS_USER $R_LIBS_USER:/usr/local/lib/R/site-library/R-lib/
ENV R_PROFILE_USER /etc/R/.Rprofile
