#!/usr/bin/env Rscript

# __author__  = "Philip Jonsson"
# __email__   = "jonssonp@mskcc.org"
# __version__ = "0.6.0"
# __status__  = "Dev"

suppressPackageStartupMessages({
    library(data.table)
    library(annotateMaf)
})

args = commandArgs(TRUE)

if (is.null(args) | length(args)<1) {
    message("Usage: filter-somatic-maf.R input.maf output.prefix")
    quit()
}

maf = args[1]
output_prefix = args[2]

add_tag = function(filter, tag) {
    ifelse(filter == 'PASS',
           tag,
           paste(filter, tag, sep = ';'))
}

maf = fread(maf, data.table = TRUE)

# Tag input MAF with filters --------------------------------------------------------------------------------------
maf[, `:=` (t_var_freq = t_alt_count/(t_alt_count+t_ref_count),
            n_var_freq = n_alt_count/(n_alt_count+n_ref_count),
            EncodeDacMapability = ifelse(is.na(EncodeDacMapability), '', EncodeDacMapability),
            RepeatMasker = ifelse(is.na(RepeatMasker), '', RepeatMasker),
            gnomAD_FILTER = ifelse(is.na(gnomAD_FILTER), 0, 1),
            FILTER = ifelse(!Custom_filters %in% c(NA, ''), Custom_filters, FILTER),
            alt_bias = t_alt_count_raw_fwd == 0 | t_alt_count_raw_rev == 0,
            ref_bias = t_depth_raw_fwd == 0 | t_depth_raw_fwd == 0
)]

maf[t_var_freq < .05, FILTER := add_tag(FILTER, 'low_vaf')]
maf[t_depth < 20, FILTER := add_tag(FILTER, 'low_t_depth')]
maf[n_depth < 10, FILTER := add_tag(FILTER, 'low_n_depth')]
maf[t_alt_count < 3, FILTER := add_tag(FILTER, 'low_t_alt_count')]
maf[n_alt_count > 3 | n_alt_count_raw > 3, FILTER := add_tag(FILTER, 'high_n_alt_count')]
maf[EncodeDacMapability != '', FILTER := add_tag(FILTER, 'mapability')]
maf[RepeatMasker != '', FILTER := add_tag(FILTER, 'repeatmasker')]
if ('non_cancer_AF_popmax' %in% names(maf)) {
    maf[non_cancer_AF_popmax > .01, FILTER := add_tag(FILTER, 'high_gnomad_pop_af')]
} else if ('AF_popmax' %in% names(maf)) {
    maf[AF_popmax > .01, FILTER := add_tag(FILTER, 'high_gnomad_pop_af')]
}
maf[PoN >= 10, FILTER := add_tag(FILTER, 'PoN')]
maf[MQ < 55 & Variant_Type != 'SNP', FILTER := add_tag(FILTER, 'low_mq')]
maf[(t_alt_count_raw > 10 & alt_bias & MuTect2 == 0) |
    (t_alt_count_raw > 10 & alt_bias & ref_bias & MQ < 40 & MuTect2 == 0), FILTER := add_tag(FILTER, 'strand_bias')]

# Filters not used:
# gnomAD_FILTER - variants considered artifacts by gnomAD's random-forest classifier

# Clean up some columns 
maf[, `:=` (Custom_filters = NULL)] 

# Tag and whitelist hotspots --------------------------------------------------------------------------------------
maf = hotspot_annotate_maf(maf)
maf = as.data.table(maf) # necessary because of the class of output from previous call
maf[Hotspot == TRUE & t_var_freq >= 0.02 & FILTER == 'low_vaf', FILTER := 'PASS'] # note: variants flagged by other filters will not be rescued by this
maf[Hotspot == TRUE & FILTER == 'low_mapping_quality', FILTER := 'PASS']
maf[Hotspot == TRUE & FILTER == 'strand_bias', FILTER := 'PASS']

# Write filtered and tagged input MAF -----------------------------------------------------------------------------
filter_maf = maf[FILTER == 'PASS']

fwrite(maf, paste0(output_prefix, '.unfiltered.maf'), sep = '\t')
fwrite(filter_maf, paste0(output_prefix, '.maf'), sep = '\t')
