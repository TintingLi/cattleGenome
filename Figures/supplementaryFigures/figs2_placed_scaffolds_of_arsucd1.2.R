################################################################################
#
# **Figure out the scaffold sequences filled by the  new Assembly.** 
#
#   1. get the unplaced sequences, and split sequences with N gaps.  
#   2. map the  unplaced seqs back to new assembly.  
#   3. check the sequence mapping results manually.
#
################################################################################


suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicAlignments))
knitr::opts_knit$set(root.dir = '~/share/litt/proj/cow/results')

setwd('~/share/litt/proj/cow/results')

#' **Final clean scaffold mapping results**  
#' 
#' + 420 scaffolds with No gaps are properly mapped to the new assembly.  
#' + 47 split scaffold sequences from 22 scaffolds with N gaps are properly mapped to the new assembly.
#' + Totally 442 scaffolds are properly mapped to the new assembly.  

scaffold.loci.final <- 
  read_tsv('./FillingGaps/Ars.ucd1.2.scaffolds/ARS_UCD1.2_scafold_loci_in_Assebmly_final.txt',
            col_names = c('seqnames', 'start', 'end', 'strand', 'qwidth', 'width', 'njuc', 
                           'scaffold', 'scaffoldLength', 'locus_chr', 'locus_chr_ctg_start', 
                          'locus_chr_ctg_end', 'locus_ctg')) %>% suppressMessages()

scaffold.loci.final$scaffold %>% duplicated() %>% table()

unique(scaffold.loci.final$scaffold) %>% length() # 467 sequences.  


grep(pattern = 'v1$', scaffold.loci.final$scaffold, value = T) %>% unique() %>%
  length() # 420 scaffolds with no N gaps.
grep(pattern = 'v1$', scaffold.loci.final$scaffold, invert = T, value = T) %>%
  unique() %>% length # 47 split regions of scaffolds containing N gaps. 

grep(pattern = 'v1$', scaffold.loci.final$scaffold, invert = T, value = T) %>% unique() %>%
  str_split(pattern = 'v1_')  %>% 
  sapply(., function(x){return(x[1])}) %>% unique() #22 scaffolds with N gaps.


#' ** Total base pairs of properly placed scaffolds**  

#' 24888609bp
filter(scaffold.loci.final, !duplicated(scaffold)) %>% 
  summarise(sum_length = sum(scaffoldLength))


scaffold.loci.final <- arrange(scaffold.loci.final[!duplicated(scaffold.loci.final$scaffold), ], desc(scaffoldLength))



#+  fig.width = 10, fig.height=14, fig.align= 'right'
par(mar= c(4,10,2,1), oma= c(4,10,2,1))
barplot(scaffold.loci.final$scaffoldLength[50:1], horiz = T, 
        names.arg = scaffold.loci.final$scaffold[50:1], las=2,
        main = 'Top 50 properly placed scaffolds according to length')













