# TAG gcr.io/dnastack-dropbox-157402/private/profyle/circos:0.69
FROM gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
	build-essential \
	libgd-dev \
	wget

RUN cpan App:cpanminus
RUN cpanm \
	Clone \
	Config::General \
	Font::TTF::Font \
	GD \
	List::MoreUtils \
	Math::Round \
	Math::Bezier \
	Math::VecStat \
	Params::Validate \
	Readonly \
	Regexp::Common \
	Set::IntSpan \
	Statistics::Basic \
	SVG \
	Text::Format

RUN wget http://circos.ca/distribution/circos-0.69-9.tgz && \
	tar -zxvf circos-0.69-9.tgz --directory /opt/ && \
	rm circos-0.69-9.tgz

RUN ln -s /opt/circos-0.69-9/bin/circos /usr/local/bin/
