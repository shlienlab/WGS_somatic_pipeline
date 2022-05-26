# TAG gcr.io/dnastack-dropbox-157402/private/profyle/gridss-purple-linx:1.3.2
FROM gridss/gridss-purple-linx:1.3.2

MAINTAINER Heather Ward <heather@dnastack.com>

# Ubuntu R GPG key is out of date
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN apt-get -qq update && \
	apt-get -qq install \
		apt-transport-https \
		ca-certificates \
		curl

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update && \
	apt-get -qq install \
		google-cloud-sdk \
		tabix \
		bcftools \
		time

ENV GRIDSS_VERSION=2.9.4

RUN rm -rf /opt/gridss
RUN mkdir -p /opt/gridss
RUN wget https://github.com/PapenfussLab/gridss/releases/download/v${GRIDSS_VERSION}/gridss-${GRIDSS_VERSION}.tar.gz && \
	tar -zxvf gridss-${GRIDSS_VERSION}.tar.gz --directory /opt/gridss && \
 	rm gridss-${GRIDSS_VERSION}.tar.gz

ENV GRIDSS_JAR=/opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar

COPY ./scripts/gridss-purple-linx.sh /opt/gridss-purple-linx/gridss-purple-linx.sh

ENTRYPOINT []
