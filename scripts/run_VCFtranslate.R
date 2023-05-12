workdir <- "/path1/studyname/"
gitdir <- "/path2/VCFtranslate/"
setwd(workdir)

# 0. prepare raw input files 
  # we need specific annotations, which can be generated using snpEff on a VCF file (+ extractFields to convert into tab-separated file)
  # additionally the column order and naming has to be exactly like this (also see example in test_data/TEST_INPUT.tsv):
  # CHROM<tab>POS<tab>REF<tab>ALT<tab>ANN[*].GENE<tab>ANN[*].EFFECT<tab>ANN[*].IMPACT<tab>ANN[*].FEATUREID<tab>ANN[*].HGVS_C<tab>ANN[*].HGVS_P<tab>patientID<tab>SOURCE
  # 19<tab>17414687<tab>A<tab>G<tab>MVB12A<tab>upstream_gene_variant<tab>MODIFIER<tab>ENST00000317040.11<tab>c.-5449A>G<tab><tab>TEST<tab>RNA

input.file <- paste0(gitdir, "/test_data/TEST_INPUT.tsv")

# variables can also be set as input parameters for this script
# args = commandArgs(trailingOnly=TRUE)
# workdir <- args[1]
# gitdir <- args[2]
# input.file <- args[3]

# 1. source and run 

source(paste0(gitdir, "/scripts/source_VCFtranslate.R"))
RunVCFtranslate(input.file)

