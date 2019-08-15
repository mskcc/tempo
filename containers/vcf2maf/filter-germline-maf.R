#!/usr/bin/env Rscript

# __author__  = "Philip Jonsson"
# __email__   = "jonssonp@mskcc.org"
# __version__ = "0.2.0"
# __status__  = "Dev"

suppressPackageStartupMessages({
    library(data.table)
    library(annotateMaf)
    library(argparse)
})

args = commandArgs(TRUE)

if (is.null(args) | length(args)<1) {
    message('Run filter-germline-maf.R --help for list of input arguments.')
    quit()
}

parser = ArgumentParser(description = 'Flag and filter somatic variants in input MAF file, output is a filtered and unfiltered but filter-tagged MAF file.')
parser$add_argument('-m', '--maf-file', required = TRUE,
                    help = 'VEP-annotated MAF file')
parser$add_argument('-o', '--output-prefix', required = TRUE,
                    help = 'Prefix of output files')
parser$add_argument('-nd', '--normal-depth', required = FALSE,
                    default = 20, help = 'Normal variant loci total depth cut-off [default %(default)s]')
parser$add_argument('-nv', '--normal-vaf', required = FALSE,
                    default = 0.35, help = 'Normal variant allele frequency cut-off [default %(default)s]')   

# Get inputs
args = parser$parse_args()
args$maf_file
output_prefix = args$output_prefix
normal_depth_cutoff = args$normal_depth
normal_vaf_cutoff = args$normal_vaf

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

maf[n_depth < normal_depth_cutoff, FILTER := add_tag(FILTER, 'low_n_depth')]
maf[!ch_gene & n_var_freq < normal_vaf_cutoff & !Variant_Classification %in% c('INS', 'DEL'), FILTER := add_tag(FILTER, 'low_n_vaf')]
maf[!ch_gene & n_var_freq < (normal_vaf_cutoff - 0.10) & Variant_Classification %in% c('INS', 'DEL'), FILTER := add_tag(FILTER, 'low_n_vaf')]
maf[ch_gene & n_var_freq < normal_vaf_cutoff & t_var_freq < .25, FILTER := add_tag(FILTER, 'ch_mutation')]
maf[t_var_freq > 3 * n_var_freq, FILTER := add_tag(FILTER, 't_in_n_contamination')]
maf[EncodeDacMapability != '', FILTER := add_tag(FILTER, 'mapability')]
maf[RepeatMasker != '', FILTER := add_tag(FILTER, 'repeatmasker')]
maf[gnomAD_FILTER == 1, FILTER := add_tag(FILTER, 'gnomad_filter')]

# Add BRCA annotation ---------------------------------------------------------------------------------------------
maf = brca_annotate_maf(maf)

# Write filtered and tagged input MAF -----------------------------------------------------------------------------
maf = as.data.table(maf) # necessary because of the class of output from previous call
filter_maf = maf[FILTER == 'PASS']

fwrite(maf, paste0(output_prefix, '.unfiltered.maf'), sep = '\t')
fwrite(filter_maf, paste0(output_prefix, '.maf'), sep = '\t')
