# VCFtranslate

## Overview
With the main goal to obtain mutated peptide sequences, mutations called from WES/WGS and RNAseq need to be introduced into the wildtype transcript DNA sequences downloaded from biomart (v92) and translated into peptide sequences. Genes are included in the analysis without exceptions regarding the transcript biotypes. For non-protein-coding transcripts, ORFs enclosing the mutation site are determined by identifying paired start and stop-codons in all 3 reading frames. The same procedure is performed for protein-coding transcripts in case of start/stop-loss/gain and frameshift mutations. Furthermore, for start/stop mutations, the coding sequence (CDS) is extended into the corresponding UTR. For mutations affecting splice donor or acceptor sites, the affected intron is included into the CDS and again checked for valid ORFs. Only mutations resulting in amino acid changes and within valid ORFs are considered. For every affected transcript, up to three ORFs enclosing the mutation site are translated into the corresponding mutated peptide sequence. Peptide sequences can then used together with the immunopeptidomics data from mass spectrometry to e.g. identify neoantigen candidates.

Input: mutations (e.g. Mutect2, Strelka) annotated with snpEFF as tab-separated text file (details below)

Output: FASTA file with mutated peptide sequences 

## Requirements
Developed and tested on R version 4.2 and 4.3.

Required packages:
- biomaRt
- data.table
- VariantAnnotation
- systemPipeR
- splitstackshape
- xml2
- Biostrings
- doParallel
- readr
- dplyr
- dbplyr
- here

## Usage
All functions are stored in scripts/source/VCFtranslate.R while the wrapper code is scripts/run_VCFtranslate.R, which needs to be opend in R or modified to be used from command line. Clone the repository and modify the input directories and file.

Input:
Specific annotations are required, which can be generated using snpEff on a VCF file (+ extractFields to convert into tab-separated file)
The column order and naming has to be exactly like this (also see example in test_data/TEST_INPUT.tsv):
```CHROM<tab>POS<tab>REF<tab>ALT<tab>ANN[*].GENE<tab>ANN[*].EFFECT<tab>ANN[*].IMPACT<tab>ANN[*].FEATUREID<tab>ANN[*].HGVS_C<tab>ANN[*].HGVS_P<tab>patientID<tab>SOURCE```
```19<tab>17414687<tab>A<tab>G<tab>MVB12A<tab>upstream_gene_variant<tab>MODIFIER<tab>ENST00000317040.11<tab>c.-5449A>G<tab><tab>TEST<tab>RNA```


