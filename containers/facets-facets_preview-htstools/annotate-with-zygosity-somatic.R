#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(binom)
  library(stringr)
})

# Usage: Rscript annotate_with_zygosity_somatic.R <input_maf> <output_maf>
args = commandArgs(trailingOnly=TRUE)

maf_file = args[1]
output_maf_file = args[2]
if ( !file.exists(maf_file) ) {
  print ("Input maf file does not exist")
  stop("Exiting")
}
 
if (length(args) < 2) {
  test_sample = "A-fb1cda9f4919-N01-IM5"
  maf_file = paste("/ifs/res/taylorlab/bandlamc/somatic_germline/somatic_germline_analysis/samples_output/", 
                   test_sample, "/", test_sample, 
                   ".germline_variants.gnomad.clinvar.facets.vep.maf", sep="")
  output_maf_file = paste("/ifs/res/taylorlab/bandlamc/somatic_germline/somatic_germline_analysis/samples_output/", 
                   test_sample, "/", test_sample, 
                   ".germline_variants.gnomad.clinvar.facets.zygosity.vep.maf", sep="")
}

t_vaf <- function(ref_copy_num, alt_copy_num, purity) { # adjusted for somatic variant
  ((purity * alt_copy_num)) / (( (purity * ref_copy_num) + (1-purity) ) + ( (purity * alt_copy_num) + (1-purity)))
}

convert_to_zygosity_string <- function(ai, loh) {
  zyg_str = ""
  if(!is.na(ai) & ai == "ALT_GAIN") {
    zygo_str = paste0("AI_", ifelse(!is.na(loh) & loh, "LOH_", ""), "ALT")
  } else if (!is.na(ai) & ai == "REF_GAIN") {
    zygo_str = paste0("AI_", ifelse(!is.na(loh) & loh, "LOH_", ""), "REF")
  } else {
    zygo_str = "None"
  }
  zygo_str
}

if (identical(maf_file, output_maf_file)){
  print ("Both input and output maf file names are identical")
  stop("Error")
}

print (paste0("Reading: ", maf_file))
maf <- suppressMessages(suppressWarnings(fread(maf_file)))

# check for required columns
req_cols <- c("t_alt_count", "t_depth", "purity", "tcn", "lcn", "cf",
              "ccf_expected_copies", "ccf_expected_copies_prob95")  # added for ccf estimate threshold -> subclonal mutations
for (col in req_cols) { if ( !(col %in% names(maf)) ) { print (paste("Missing required column in maf: ", col, sep=""))} }

# initialize zygosity related variables
# 
maf <-
  maf %>%
  mutate(num_ref_copies = NA,
         num_alt_copies = NA,
         expected_t_alt_freq_lower = NA,
         expected_t_alt_freq_upper = NA,
         tumor_vaf_cn_concordance = NA,
         expected_t_alt_freq_lower_99 = NA,
         expected_t_alt_freq_upper_99 = NA,
         tumor_vaf_cn_concordance_99 = NA,
         allelic_imbalance = NA,
         loss_of_heterozygosity = NA, 
         zygosity_flag = "None")

# calculate zygosity for each variant in the maf
for (idx in 1:dim(maf)[1]){

  ## This optional column indicates whether facets has been run for this sample
  # if (is.na(maf$has_FACETS[idx]) | maf$has_FACETS[idx] == 0){
  #   next
  # }

  lcn = maf$lcn[idx]
  tcn = maf$tcn[idx]
  purity = maf$purity[idx]
  t_depth = maf$t_depth[idx]
  t_alt_count = maf$t_alt_count[idx]
  t_alt_freq = ifelse('t_alt_freq' %in% names(maf), maf$t_alt_freq[idx], maf$t_alt_count[idx]/maf$t_depth[idx])
  # n_depth = maf$n_depth[idx]
  # n_alt_freq = maf$n_alt_freq[idx]

  # if ( is.na(purity) | is.na(t_alt_freq) | is.infinite(purity) | n_alt_freq > 0.75 | n_depth < 50 ) {
  if ( is.na(purity) | is.na(t_alt_freq) | is.infinite(purity) | (is.na(tcn) & is.na(lcn))) {
    #skipping these, we cant do anything
      maf$zygosity_flag[idx] = 'Indeterminate'
    if (is.infinite(purity)) {
      maf$Purity[idx] = NA
    }
    next
  }

  mcn = tcn - lcn

  if (maf$ccf_expected_copies_upper[idx] < 0.80 & maf$ccf_expected_copies_prob95[idx] < 0.80) {
    maf$zygosity_flag[idx] = 'subclonal'
    next
  }

  # shouldn't happen; a legacy check 
  if ( t_alt_count > t_depth ) { 
    print ("WARNING: Doesn't make sense but we see this with the DMP pipeline.  Shouldn't happen anymore because we are doing maf-fill.") 
  }
  
  # reference alignment bias correction factor
  ref_aln_bias_factor = 1 # not relevant for somatic variants, set to 1

  # allelic imbalance (values: REF-GAIN/ALT-GAIN)
  # flag as AI if the point estimate for t_alt_freq is above the lower CI of the 2-1 or 1-2 CN states 
  # for somatic variants, "imbalance treshold" is the heterozygous state
  if (t_alt_freq >= purity/2) {
    imbalanced_vaf_ci = binom.confint( min(t_depth, t_vaf(1, 2, purity) * t_depth * ref_aln_bias_factor), 
                                       t_depth, methods="wilson")
    if (imbalanced_vaf_ci$lower <= t_alt_freq){
      maf$allelic_imbalance[idx] = "ALT_GAIN"
    }
  } else {
    imbalanced_vaf_ci = binom.confint( min(t_depth, t_vaf(2, 1, purity) * t_depth * ref_aln_bias_factor), 
                                       t_depth, methods="wilson")
    if (imbalanced_vaf_ci$upper >= t_alt_freq){
      maf$allelic_imbalance[idx] = "REF_GAIN"
    }
  }
  
  if (is.na(maf$allelic_imbalance[idx])) { # tcn == 0 and no imbalance
      next
  }  else if (is.na(lcn)) {
      if (maf$allelic_imbalance[idx] == "ALT_GAIN") { # is.na(lcn) and imbalance favoring ALT
          expected_vaf_loh = t_vaf(0, tcn, purity) * ref_aln_bias_factor
          expected_vaf_noloh = t_vaf(1, tcn-1, purity) * ref_aln_bias_factor
          thresh_vaf_loh = (binom.confint(min(t_depth, t_vaf(0, tcn, purity) * t_depth * ref_aln_bias_factor), 
                                          t_depth, 
                                          methods="wilson"))$lower
          if (abs(expected_vaf_loh - t_alt_freq) < abs(expected_vaf_noloh - t_alt_freq) &
              t_alt_freq > thresh_vaf_loh) {
              ref_cn = 0
              alt_cn = tcn
              maf$loss_of_heterozygosity[idx] = TRUE
              maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
              next
          } else {
              maf$loss_of_heterozygosity[idx] = FALSE
              maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
              next
          }
      } else if (maf$allelic_imbalance[idx] == "REF_GAIN") { # is.na(lcn) and imbalance favoring REF
          expected_vaf_loh = t_vaf(tcn, 0, purity) * ref_aln_bias_factor
          expected_vaf_noloh = t_vaf(tcn-1, 1, purity) * ref_aln_bias_factor   
          thresh_vaf_loh = (binom.confint(min(t_depth, t_vaf(tcn, 0, purity) * t_depth * ref_aln_bias_factor), 
                                          t_depth, 
                                          methods="wilson"))$upper
          if (abs(expected_vaf_loh - t_alt_freq) < abs(expected_vaf_noloh - t_alt_freq) &
              t_alt_freq < thresh_vaf_loh) {
              ref_cn = tcn
              alt_cn = 0
              maf$loss_of_heterozygosity[idx] = TRUE
              maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
              next
          } else {
              maf$loss_of_heterozygosity[idx] = FALSE
              maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
              next
          }
      } else { next } # is.na(lcn) and no imbalance
  } else if (tcn == 0) {
      mcn = tcn = 1 # assumption that a copy exists
      if (maf$allelic_imbalance[idx] == "ALT_GAIN") { # tcn == 0 and imbalance favoring ALT
          expected_vaf_loh = t_vaf(0, tcn, purity) * ref_aln_bias_factor
          expected_vaf_noloh = t_vaf(1, tcn-1, purity) * ref_aln_bias_factor   
          thresh_vaf_loh = (binom.confint(min(t_depth, t_vaf(0, tcn, purity) * t_depth * ref_aln_bias_factor), 
                                          t_depth, 
                                          methods="wilson"))$lower
          if (abs(expected_vaf_loh - t_alt_freq) < abs(expected_vaf_noloh - t_alt_freq) &
              t_alt_freq > thresh_vaf_loh) {
              ref_cn = 0
              alt_cn = tcn
              maf$loss_of_heterozygosity[idx] = TRUE
              maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
              next
          } else {
              maf$loss_of_heterozygosity[idx] = FALSE
              maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
              next
          }
      } else if (maf$allelic_imbalance[idx] == "REF_GAIN") { # tcn == 0 and imbalance favoring REF
          expected_vaf_loh = t_vaf(tcn, 0, purity) * ref_aln_bias_factor
          expected_vaf_noloh = t_vaf(tcn-1, 1, purity) * ref_aln_bias_factor   
          thresh_vaf_loh = (binom.confint(min(t_depth, t_vaf(tcn, 0, purity) * t_depth * ref_aln_bias_factor), 
                                          t_depth, 
                                          methods="wilson"))$upper
          if (abs(expected_vaf_loh - t_alt_freq) < abs(expected_vaf_noloh - t_alt_freq) &
              t_alt_freq < thresh_vaf_loh) {
              ref_cn = tcn
              alt_cn = 0
              maf$loss_of_heterozygosity[idx] = TRUE
              maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
              next
          } else {
              maf$loss_of_heterozygosity[idx] = FALSE
              maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
              next
          }
      } 
  } 
  
  # Check if any of the required columns are NA
  # expected VAF of the variant given that the variant copy number in tumor == mcn
  expected_vaf_mcn = t_vaf(lcn, mcn, purity) * ref_aln_bias_factor  
  
  # expected VAF of the variant given that the variant copy number in tumor == lcn
  expected_vaf_lcn = t_vaf(mcn, lcn, purity) * ref_aln_bias_factor  
  
  # Infer the number of mutant copies (this is by no means a precise estimate, and, should be interpreted with caution):
  #  -- estimate copy number by determining whether t_alt_freq is closer to expected_vaf_mcn or expected_vaf_lcn
  if (abs(expected_vaf_mcn - t_alt_freq) < abs(expected_vaf_lcn - t_alt_freq)) {
    expected_vaf = expected_vaf_mcn 
    ref_cn = lcn
    alt_cn = mcn
  } else {
    expected_vaf = expected_vaf_lcn
    ref_cn = mcn
    alt_cn = lcn
  }
  maf$num_ref_copies[idx] = ref_cn
  maf$num_alt_copies[idx] = alt_cn
  
  expected_vaf_ci = binom.confint( min(t_depth, expected_vaf * t_depth * ref_aln_bias_factor), 
                                   t_depth, methods="wilson" )
  
  # TRUE/FALSE.  Concordant if the the observed tumor_vaf given the purity and num_alt_copies falls within 
  # the expected range determined using binomial confidence intervals.  
  if (t_alt_freq >= expected_vaf_ci$lower && t_alt_freq <= expected_vaf_ci$upper ) {
    maf$tumor_vaf_cn_concordance[idx] = TRUE
  } else {
    maf$tumor_vaf_cn_concordance[idx] = FALSE
  }
  
  maf$expected_t_alt_freq_lower[idx] = round(expected_vaf_ci$lower, 3)
  maf$expected_t_alt_freq_upper[idx] = round(expected_vaf_ci$upper, 3)
  
  expected_vaf_ci_99 = binom.confint( min(t_depth, expected_vaf * t_depth * ref_aln_bias_factor), 
                                   t_depth, conf.level = 0.99, methods="wilson" )
  
  # TRUE/FALSE.  Concordant if the the observed tumor_vaf given the purity and num_alt_copies falls within 
  # the expected range determined using binomial confidence intervals.  
  if ( grepl("GAIN", maf$allelic_imbalance[idx]) ) {
      # Call as "LOH" if ref_cn/alt_cn is 0 or if the tumor variant allele frequencies are concordant with
      # a complete loss of ref/alt copies
      if ((grepl("ALT_GAIN", maf$allelic_imbalance[idx]) && ref_cn == 0) | 
          (grepl("REF_GAIN", maf$allelic_imbalance[idx]) && alt_cn == 0)) {
          maf$loss_of_heterozygosity[idx] = TRUE
      } else if ( maf$allelic_imbalance[idx] == "ALT_GAIN" ) {
          thresh_vaf = (binom.confint( min(t_depth, t_vaf(0, alt_cn, purity) * t_depth * ref_aln_bias_factor), 
                         t_depth, 
                         methods="wilson"))$lower
          if ( t_alt_freq >= thresh_vaf) {
              maf$loss_of_heterozygosity[idx] = TRUE
          }
      } else {
          thresh_vaf = (binom.confint( min(t_depth, t_vaf(ref_cn, 0, purity) * t_depth * ref_aln_bias_factor), 
                                            t_depth, 
                                            methods="wilson"))$upper
          if ( t_alt_freq <= thresh_vaf) {
              maf$loss_of_heterozygosity[idx] = TRUE
          }          
      }
  }
  maf$zygosity_flag[idx] = convert_to_zygosity_string(maf$allelic_imbalance[idx], maf$loss_of_heterozygosity[idx])
}

write.table(maf, file=output_maf_file, sep="\t", row.names=F, quote=F)
