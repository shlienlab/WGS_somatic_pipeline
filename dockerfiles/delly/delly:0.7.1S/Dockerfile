# TAG gcr.io/dnastack-dropbox-157402/private/profyle/delly:0.7.1S
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

RUN apt-get -qq update && \
	apt-get -qq install \
		libhts-dev

# Shlien lab's version of Delly and boost libs
COPY ./modular-boost.tar.gz /opt/
RUN tar -zxvf /opt/modular-boost.tar.gz --directory /opt/ && \
	rm /opt/modular-boost.tar.gz

COPY ./scripts/delly_v0.7.1S_CentOS5.4_x86_64bit /usr/local/bin/
RUN ln -s /usr/local/bin/delly_v0.7.1S_CentOS5.4_x86_64bit /usr/local/bin/delly


ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/modular-boost/stage/lib
