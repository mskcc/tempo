library(argparse)

parser=ArgumentParser()
parser$add_argument("-f", "--facets_maf",type="character",help="facets input maf file")
parser$add_argument("-n", "--neoantigen_maf",type="character",help="neoantigen input maf file")
parser$add_argument("-o", "--output_file",type="character",default="merged.output.maf",help="output maf file name")

args=parser$parse_args()

FACETS_MAF=args$facets_maf
NEO_MAF=args$neoantigen_maf
OUTPUT_FILE=args$output_file

library(data.table)
facets_maf=fread(FACETS_MAF)
neo_maf=fread(NEO_MAF)

neoantigen_cols = c("neo_maf_identifier_key", "neo_best_icore_peptide", "neo_best_rank", "neo_best_binding_affinity", "neo_best_binder_class", "neo_best_is_in_wt_peptidome", "neo_best_algorithm", "neo_best_hla_allele","neo_n_peptides_evaluated", "neo_n_strong_binders", "neo_n_weak_binders")       

final_maf=cbind(facets_maf, neo_maf[, ..neoantigen_cols])
fwrite(final_maf, OUTPUT_FILE, sep="\t")
