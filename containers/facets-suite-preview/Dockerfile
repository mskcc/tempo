FROM rocker/tidyverse:3.6.1

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
	contributor="Yixiao Gong (gongy@mskcc.org)" \
	contributor="Nikhil Kumar (kumarn1@mskcc.org)" \
	contributor="Philip Jonsson (jonssonp@mskcc.org)" \
	version.image="0.0.1-rebuild" \
	version.facets_suite="2.0.8" \
	version.facets="0.5.14" \
	version.alpine="3.8" \
	version.pctGCdata="0.2.0" \
	source.facets="https://github.com/mskcc/facets/archive/v0.5.14.tar.gz" 

ENV FACETS_SUITE_VERSION 2.0.8
ENV FACETS_VERSION 0.5.14
ENV FACETS_PREVIEW_VERSION 2.1.4
ENV PCTGCDATA 0.2.0

# Requirements
RUN apt-get update \
	&& apt-get install -y \
	g++ \
	tar \
	bzip2 \
	libbz2-dev \
	libc6-dev \
	libxt-dev \
	liblzma-dev \
	libgtk2.0-dev \
	libcairo2-dev \
	xvfb \
	xauth \
	xfonts-base
RUN apt-get install -y xdg-utils --fix-missing
RUN apt-get install -y inotify-tools=3.14-2

RUN R -e "install.packages(c('Cairo','argparse','gridExtra', 'binom', 'BiocManager', 'diptest', 'egg','shinyWidgets', 'shinyjs', 'rhandsontable', 'doParallel', 'configr', 'R.utils'), repos='http://cran.us.r-project.org')" \
	&& R -e "BiocManager::install('rtracklayer')"

# Install FACETS, pctGCdata and facets-suite
RUN cd /tmp \
	&& wget https://github.com/mskcc/facets-suite/archive/${FACETS_SUITE_VERSION}.tar.gz -O facets-suite-${FACETS_SUITE_VERSION}.tar.gz \
	&& wget https://github.com/mskcc/facets/archive/v${FACETS_VERSION}.tar.gz -O facets-v${FACETS_VERSION}.tar.gz \
	&& wget https://github.com/mskcc/pctGCdata/archive/v${PCTGCDATA}.tar.gz \
	&& git clone --single-branch --branch master https://github.com/taylor-lab/facets-preview.git facets-preview-master \
	&& tar xvzf facets-v${FACETS_VERSION}.tar.gz \
	&& tar xvzf v${PCTGCDATA}.tar.gz \
	&& tar xvzf facets-suite-${FACETS_SUITE_VERSION}.tar.gz \
	&& cd /tmp/pctGCdata-${PCTGCDATA} \
	&& R CMD INSTALL . \
	&& cd /tmp/facets-${FACETS_VERSION} \
	&& R CMD INSTALL . \
	&& cd /tmp/facets-suite-${FACETS_SUITE_VERSION} \
	&& R CMD INSTALL . \
	# correct shebang line
	&& sed -i "s/opt\/common\/CentOS_6-dev\/R\/R-3.2.2\//usr\//g" *.R \
	# copy execs to /usr/bin/facets-suite
	&& mkdir -p /usr/bin/facets-suite/ \
	&& cp -r /tmp/facets-suite-${FACETS_SUITE_VERSION}/* /usr/bin/facets-suite/ \
	&& cd /tmp/facets-preview-master \
	&& R CMD INSTALL . \
	&& mkdir -p /usr/bin/facets-preview/ \
	&& cp -r /tmp/facets-preview-master/* /usr/bin/facets-preview/ \
	# clean up
	&& rm -rf /var/cache/apk/* /tmp/* # update to using release version when 2.1.5+ is released.

ENV PYTHONNOUSERSITE set
ENV FACETS_OVERRIDE_EXITCODE set

COPY impact_juno_config.json /usr/bin/facets-preview/
RUN chmod 744 /usr/bin/facets-preview/impact_juno_config.json
EXPOSE 3838

