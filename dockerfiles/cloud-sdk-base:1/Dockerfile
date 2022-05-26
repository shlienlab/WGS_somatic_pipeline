# TAG gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1
FROM ubuntu:xenial

RUN apt-get -qq update && \
	apt-get -qq install \
		apt-transport-https \
		ca-certificates \
		gnupg \
		curl

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

RUN apt-get -qq update && \
	apt-get -qq install \
		google-cloud-sdk
