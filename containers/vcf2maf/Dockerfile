from cmopipeline/vcf2maf:vep88_1.2.7

LABEL maintainer="Christopher Allan Bolipata (bolipatc@mskcc.org)" \
      contributor="Nikhil Kumar (kumarn1@mskcc.org)" \
      contributor="Philip Jonsson (jonssonp@mskcc.org)" \
      contributor="Anne Marie Noronha (noronhaa@mskcc.org)" \
      version.image="1.3.0" \
      version.vcf2maf="1.6.17" \
      version.vep="88" \
      version.htslib="1.9" \
      version.bcftools="1.9" \
      version.samtools="1.9" \
      version.perl="5.26.2-r1" \
      version.alpine="3.8" \
      version.oncokb_annotator="1.1.0" \
      version.filter_somatic_maf="0.7.0" \
      version.filter_germline_maf="0.2.2" \
      version.annotateMaf="1.0.3" \
      source.vcf2maf="https://github.com/mskcc/vcf2maf/releases/tag/v1.6.17" \
      source.vep="http://dec2016.archive.ensembl.org/info/docs/tools/vep/script/vep_download.html#versions" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.9" \
      source.bcftools="https://github.com/samtools/bcftools/releases/tag/1.9" \
      source.samtools="https://github.com/samtools/samtools/releases/tag/1.9" \
      source.oncokb_annotator="https://github.com/oncokb/oncokb-annotator/releases/tag/v1.1.0" \
      source.annotateMaf="https://github.com/taylor-lab/annotateMaf"

ENV ANNOTATE_MAF 1.0.3

RUN R -e "devtools::install_github('taylor-lab/annotateMaf', ref = '${ANNOTATE_MAF}', upgrade = 'never')"
   
# Add filter script
RUN rm /usr/bin/filter-somatic-maf.R
RUN rm /usr/bin/filter-germline-maf.R

COPY filter-somatic-maf.R /usr/bin
COPY filter-germline-maf.R /usr/bin
RUN chmod 755 /usr/bin/filter-somatic-maf.R
RUN chmod 755 /usr/bin/filter-germline-maf.R

