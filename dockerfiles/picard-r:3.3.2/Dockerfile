# TAG gcr.io/dnastack-dropbox-157402/private/profyle/picard-r:3.3.2
FROM gcr.io/dnastack-dropbox-157402/private/profyle/picard:1.108

MAINTAINER Heather Ward <heather@dnastack.com>

RUN wget https://cran.r-project.org/bin/linux/ubuntu/xenial/r-base-core_3.3.2-1xenial0_amd64.deb && \
	apt-get -qq install ./r-base-core_3.3.2-1xenial0_amd64.deb

RUN wget https://cran.r-project.org/bin/linux/ubuntu/xenial/r-recommended_3.3.2-1xenial0_all.deb && \
	apt-get -qq install ./r-recommended_3.3.2-1xenial0_all.deb

RUN wget https://cran.r-project.org/bin/linux/ubuntu/xenial/r-base_3.3.2-1xenial0_all.deb && \
	apt-get -qq install ./r-base_3.3.2-1xenial0_all.deb
