#!/usr/bin/env Rscript

# __author__  = "Philip Jonsson"
# __email__   = "jonssonp@mskcc.org"
# __version__ = "0.6.0"
# __status__  = "Dev"

suppressPackageStartupMessages({
    library(data.table)
    library(annotateMaf)
    library(argparse)
})

args = commandArgs(TRUE)

if (is.null(args) | length(args)<1) {
    message('Run filter-somatic-maf.R --help for list of input arguments.')
    quit()
}

parser = ArgumentParser(description = 'Flag and filter somatic variants in input MAF file, output is a filtered and unfiltered but filter-tagged MAF file.')
parser$add_argument('-m', '--maf-file', required = TRUE,
                    help = 'VEP-annotated MAF file')
parser$add_argument('-o', '--output-prefix', required = TRUE,
                    help = 'Prefix of output files')
parser$add_argument('-tv', '--tumor-vaf', required = FALSE,
                    default = 0.05, help = 'Tumor variant allele frequency cut-off [default %(default)s]')
parser$add_argument('-td', '--tumor-depth', required = FALSE,
                    default = 20, help = 'Tumor variant loci total depth cut-off [default %(default)s]')
parser$add_argument('-tc', '--tumor-count', required = FALSE,
                    default = 3, help = 'Tumor variant read support cut-off [default %(default)s]')                    
parser$add_argument('-nd', '--normal-depth', required = FALSE,
                    default = 10, help = 'Normal variant loci total depth cut-off [default %(default)s]')
parser$add_argument('-nc', '--normal-count', required = FALSE,
                    default = 3, help = 'Normal variant read support cut-off [default %(default)s]')
parser$add_argument('-gaf', '--gnomad-allele-frequency', required = FALSE,
                    default = 0.01, help = 'gnomAD allele frequency cut-off [default %(default)s]')
parser$add_argument('-pon', '--normal-panel-count', required = FALSE,
                    default = 10, help = 'Panel of normals count cut-off [default %(default)s]')
                
# Get inputs
args = parser$parse_args()
maf = args$maf_file
output_prefix = args$output_prefix
tumor_vaf_cutoff = args$tumor_vaf
tumor_depth_cutoff = args$tumor_depth
tumor_readcount_cutoff = args$tumor_count
normal_depth_cutoff = args$normal_depth
normal_readcount_cutoff = args$normal_count
gnomad_af_cutoff = args$gnomad_allele_frequency
pon_cutoff = args$normal_panel_count

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
            alt_bias = t_depth_raw > 5 & (t_alt_count_raw_fwd == 0 | t_alt_count_raw_rev == 0),
            ref_bias = t_depth_raw > 5 & (t_depth_raw_fwd == 0 | t_depth_raw_fwd == 0)
)]

# Apply filters to FILTER column
maf[t_var_freq < tumor_vaf_cutoff, FILTER := add_tag(FILTER, 'low_vaf')]
maf[t_depth < tumor_depth_cutoff, FILTER := add_tag(FILTER, 'low_t_depth')]
maf[t_alt_count < tumor_readcount_cutoff, FILTER := add_tag(FILTER, 'low_t_alt_count')]
maf[n_depth < normal_depth_cutoff, FILTER := add_tag(FILTER, 'low_n_depth')]
maf[n_alt_count > normal_readcount_cutoff | n_alt_count_raw > normal_readcount_cutoff, FILTER := add_tag(FILTER, 'high_n_alt_count')]
maf[EncodeDacMapability != '', FILTER := add_tag(FILTER, 'mapability')]
maf[RepeatMasker != '', FILTER := add_tag(FILTER, 'repeatmasker')]
if ('non_cancer_AF_popmax' %in% names(maf)) {
    maf[non_cancer_AF_popmax > gnomad_af_cutoff, FILTER := add_tag(FILTER, 'high_gnomad_pop_af')]
} else if ('AF_popmax' %in% names(maf)) {
    maf[AF_popmax > gnomad_af_cutoff, FILTER := add_tag(FILTER, 'high_gnomad_pop_af')]
}
maf[PoN >= pon_cutoff, FILTER := add_tag(FILTER, 'PoN')]
maf[MQ < 55 & Variant_Type != 'SNP', FILTER := add_tag(FILTER, 'low_mapping_quality')]
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
