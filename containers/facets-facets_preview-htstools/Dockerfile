FROM cmopipeline/facets-htstools:3.0.3

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
      version.image="0.0.1" \
      version.facets-htstools="3.0.3" \
      version.facets_preview="2.1.4" 
      
ENV FACETS_PREVIEW_VERSION 2.1.4

RUN R -e "install.packages(c('shinyWidgets', 'shinyjs', 'rhandsontable', 'doParallel', 'configr', 'R.utils'), repos='http://cran.us.r-project.org')" 

RUN cd /tmp \
    && wget https://github.com/taylor-lab/facets-preview/archive/${FACETS_PREVIEW_VERSION}.tar.gz -O facets-preview-${FACETS_PREVIEW_VERSION}.tar.gz \
    && tar xvzf facets-preview-${FACETS_PREVIEW_VERSION}.tar.gz \
    && cd facets-preview-${FACETS_PREVIEW_VERSION} \
    && R CMD INSTALL . \
    && mkdir -p /usr/bin/facets-preview/ \
    && cp -r /tmp/facets-preview-${FACETS_PREVIEW_VERSION}/* /usr/bin/facets-preview/ \
    && rm -rf /tmp/*

COPY tempo_config.json /usr/bin/facets-preview/
RUN chmod 744 /usr/bin/facets-preview/tempo_config.json

