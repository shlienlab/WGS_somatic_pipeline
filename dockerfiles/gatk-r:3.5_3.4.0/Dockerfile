# TAG gcr.io/dnastack-dropbox-157402/private/profyle/gatk-r:3.5_3.4.0
FROM gcr.io/dnastack-dropbox-157402/private/profyle/bcftools-r:3.4.0

RUN apt-get -qq update && \
	apt-get -qq install \
	software-properties-common \
	bzip2 \
	tabix

RUN add-apt-repository ppa:openjdk-r/ppa
RUN apt-get -qq update && \
	apt-get -qq install \
	openjdk-7-jdk

RUN update-alternatives --set java /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java

RUN gsutil cp gs://gatk-software/package-archive/gatk/GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2 . && \
	tar -jxvf GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2 --directory /opt/ && \
	rm GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2 && \
	rm -rf /opt/resources

ENV GATK /opt/GenomeAnalysisTK.jar
