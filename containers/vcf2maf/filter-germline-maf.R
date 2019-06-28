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
    message("Usage: filter-germline-maf.R input.maf")
    quit()
}

maf = args[1]
output_prefix = args[2]

add_tag = function(filter, tag) {
    ifelse(filter == 'PASS',
           tag,
           paste(filter, tag, sep = ';'))
}

ch_genes = c("ASXL1", "ATM", "BCOR", "CALR", "CBL", "CEBPA", "CREBBP", "DNMT3A", "ETV6", "EZH2", "FLT3", "GNAS",
             "IDH1", "IDH2", "JAK2", "KIT", "KRAS", "MPL", "MYD88", "NF1", "NPM1", "NRAS", "PPM1D", "RAD21", "RUNX1", "SETD2", "SF3B1", "SH2B3", "SRSF2", "STAG2", "STAT3", "TET2", "TP53", "U2AF1", "WT1", "ZRSR2")

maf = fread(maf, data.table = TRUE)

# Tag input MAF with filters --------------------------------------------------------------------------------------
maf[, `:=` (t_var_freq = t_alt_count/(t_alt_count+t_ref_count),
            n_var_freq = n_alt_count/(n_alt_count+n_ref_count),
            EncodeDacMapability = ifelse(is.na(EncodeDacMapability), '', EncodeDacMapability),
            RepeatMasker = ifelse(is.na(RepeatMasker), '', RepeatMasker),
            gnomAD_FILTER = ifelse(is.na(gnomAD_FILTER), 0, 1),
            ch_gene = Hugo_Symbol %in% ch_genes
)]

maf[n_depth < 20, FILTER := add_tag(FILTER, 'low_n_depth')]

maf[!ch_gene & n_var_freq < .35 & !Variant_Classification %in% c('INS', 'DEL'), FILTER := add_tag(FILTER, 'low_n_vaf')]
maf[!ch_gene & n_var_freq < .25 & Variant_Classification %in% c('INS', 'DEL'), FILTER := add_tag(FILTER, 'low_n_vaf')]
maf[ch_gene & n_var_freq < .35 & t_var_freq < .25, FILTER := add_tag(FILTER, 'ch_mutation')]
maf[t_var_freq > 3 * n_var_freq, FILTER := add_tag(FILTER, 't_in_n_contamination')]
maf[EncodeDacMapability != '', FILTER := add_tag(FILTER, 'mapability')]
maf[RepeatMasker != '', FILTER := add_tag(FILTER, 'repeatmasker')]
maf[gnomAD_FILTER == 1, FILTER := add_tag(FILTER, 'gnomad_filter')]

# Filters not used:
# gnomAD_FILTER - variants considered artifacts by gnomAD's random-forest classifier

# Add BRCA annotation ---------------------------------------------------------------------------------------------
maf = brca_annotate_maf(maf)

# Write filtered and tagged input MAF -----------------------------------------------------------------------------
maf = as.data.table(maf) # necessary because of the class of output from previous call
filter_maf = maf[FILTER == 'PASS']

fwrite(maf, paste0(output_prefix, '.unfiltered.maf'), sep = '\t')
fwrite(filter_maf, paste0(output_prefix, '.maf'), sep = '\t')
