#!/usr/bin/env Rscript

# __author__  = "Philip Jonsson"
# __email__   = "jonssonp@mskcc.org"
# __version__ = "0.1.0"
# __status__  = "Dev"

suppressPackageStartupMessages({
    library(data.table)
})

args = commandArgs(TRUE)

if (is.null(args) | length(args)<1) {
    message("Usage: filter-maf.R input.maf")
    quit()
}

maf = args[1]
output1 = maf
output2 = gsub('.unfiltered.maf$', '.maf', maf)

add_tag = function(filter, tag) {
    ifelse(filter == 'PASS',
           tag,
           paste(filter, tag, sep = ';'))
}

# Tag input MAF with filters --------------------------------------------------------------------------------------
maf = fread(maf)

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
maf[n_alt_count > 3, FILTER := add_tag(FILTER, 'high_n_alt_count')]
maf[EncodeDacMapability != '', FILTER := add_tag(FILTER, 'mapability')]
maf[RepeatMasker != '', FILTER := add_tag(FILTER, 'repeatmasker')]
maf[non_cancer_AF_popmax > .01, FILTER := add_tag(FILTER, 'high_gnomad_pop_af')]
maf[PoN >= 10, FILTER := add_tag(FILTER, 'PoN')]

# Filters not used:
# gnomAD_FILTER - variants considered artifacts by gnomAD's random-forest classifier

filter_maf = maf[FILTER == 'PASS']

# Write filtered and tagged input MAF -----------------------------------------------------------------------------
fwrite(maf, output1, sep = '\t')
fwrite(filter_maf, output2, sep = '\t')
