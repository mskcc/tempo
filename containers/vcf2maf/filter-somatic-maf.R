#!/usr/bin/env Rscript

# __author__  = "Philip Jonsson"
# __email__   = "jonssonp@mskcc.org"
# __version__ = "0.4.0"
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
            gnomAD_FILTER = ifelse(is.na(gnomAD_FILTER), 0, 1)
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

# Filters not used:
# gnomAD_FILTER - variants considered artifacts by gnomAD's random-forest classifier

# Tag hotspots ----------------------------------------------------------------------------------------------------
maf = hotspot_annotate_maf(maf)

# Write filtered and tagged input MAF -----------------------------------------------------------------------------
maf = as.data.table(maf) # necessary because of the class of output from previous call
filter_maf = maf[FILTER == 'PASS']

fwrite(maf, paste0(output_prefix, '.unfiltered.maf'), sep = '\t')
fwrite(filter_maf, paste0(output_prefix, '.maf'), sep = '\t')
