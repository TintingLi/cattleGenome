#' ---
#' title: "Fetch V nt sequences of the annotated immuno genomic loci"
#' author: "LiTingting(ting67@126.com)"
#' output:
#'   html_document:
#'     number_sections: true
#'     highlight: pygments
#'     theme: cerulean
#'     toc: yes
#' ---



suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(Biostrings))

setwd('~/Desktop/projects/cow/results/')
knitr::opts_knit$set(root.dir = '~/Desktop/projects/cow/results/')




#' ## IG immuno loci
#' **Note that for IG loci, genes were mapped directly to the fetched sequences.**  
#' 
#' ### fetch IGHV sequence of the new assembly (rev in genome)  
igh <- readDNAStringSet('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/IGH/mapToIGH_revComp/bt2idx/IGH_ctg114_67750912_68414741_revcomp.fa')
igh.anno <- read_csv('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/IGH/mapToIGH_revComp/IGH_manuallyChecked.csv') %>% suppressMessages()
igh.v <- igh.anno[grepl('IGHV', igh.anno$annotation),  ]
ighv.fa <- Views(igh[[1]], start = igh.v$start, width = igh.v$length)
names(ighv.fa) <- paste0(igh.v$annotation, '_', igh.v$functionality, '_ncba.01')
ighv.fa <- DNAStringSet(ighv.fa)
ighv.fa
writeXStringSet(ighv.fa, './phylogen/fa/ncba.01/bovine_cattle.ncba01.IGHV.fa')


#' ### fetch IGLV sequence of the new assembly (fwd in genome)  

igl <- readDNAStringSet('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/IGL/mapToIGL_locus/bt2idx/IGL_chr17_72_72.7.fa')
igl.anno <- read_csv('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/IGL/mapToIGL_locus/xiatian/IGL_mapResults-20220130_cleanup.csv') %>% suppressMessages()
igl.v <- igl.anno[grepl('IGLV', igl.anno$annotation),  ]
iglv.fa <- Views(igl[[1]], start = igl.v$start, width = igl.v$length)
names(iglv.fa) <- paste0(igl.v$annotation, '_', igl.v$functionality, '_ncba.01')
iglv.fa <- DNAStringSet(iglv.fa)
iglv.fa
writeXStringSet(iglv.fa, './phylogen/fa/ncba.01/bovine_cattle.ncba01.IGLV.fa')



#' ###  fetch IGKV sequence of the new assembly (fwd in genome)  

IGK <- readDNAStringSet('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/IGK/mapToIGK_locus/bt2idx/IGK_chr11_47.1_47.5.fa')
IGK.anno <- read_csv('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/IGK/mapToIGK_locus/IGK_mapResults_clean.csv') %>% suppressMessages()

IGK.v <- IGK.anno[grepl('IGKV', IGK.anno$annotation),  ]
IGKv.fa <- Views(IGK[[1]], start = IGK.v$start, width = IGK.v$length)
names(IGKv.fa) <- paste0(IGK.v$annotation, '_', IGK.v$functionality, '_ncba.01')

IGKv.fa <- DNAStringSet(IGKv.fa)
IGKv.fa
writeXStringSet(IGKv.fa, './phylogen/fa/ncba.01/bovine_cattle.ncba01.IGKV.fa')



#' ## TR immuno loci  
#' 
#' ** Note that for TR loci, genes were mapped to the new assembly.**  
#' **So, if the locus was rev in the genome, the gene sequences should be reversecomplemented after fetch.**
#' 
#' ###  fetch TRA_D sequence of the new assembly (rev in genome)

TRA_D <- readDNAStringSet('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/TRA_TRD/TRA_D_chr10_22.0_26.0Mb.fa')
trad.anno <- read_csv('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/TRA_TRD/TRA_TRD_mapResults_clean_merged.csv')%>% suppressMessages()

trad.v.anno <- trad.anno[grepl('TRAV|TRDV', trad.anno$annotation),  ]
trad.v <- Views(TRA_D[[1]], start = (trad.v.anno$start - 22000000 +1), width = trad.v.anno$length) %>% reverseComplement()
names(trad.v) <- paste0(trad.v.anno$annotation, '_', trad.v.anno$functionality, '_ncba.01')


trad.v <- DNAStringSet(trad.v)
#' **TRDV3 is in different direction**  
trad.v$TRDV3_F <- reverseComplement(trad.v$TRDV3_F)
trad.v
writeXStringSet(DNAStringSet(trad.v), './phylogen/fa/ncba.01/bovine_cattle.ncba01.TRA_TRD_V.fa')



#' ###  fetch TRB sequence of the new assembly (FWD in genome)

TRB <- readDNAStringSet('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/TRB/TRB_chr4_105.4_106.3Mb.fa')
TRB.anno <- read_csv('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/TRB/xiatian/TRB_GENEDB_clean.csv')%>% suppressMessages()


TRB.anno$length <-  str_split(TRB.anno$gene, '_') %>% sapply( function(x){return(x[4])}) %>% sub('nt$', '', x=.) %>% as.numeric()


TRB.v.anno <- TRB.anno[grepl('TRBV', TRB.anno$annotation),  ]
TRB.v <- Views(TRB[[1]], start = (TRB.v.anno$start - 105400000 +1), width = TRB.v.anno$length)

names(TRB.v) <- paste0(TRB.v.anno$annotation, '_', TRB.v.anno$functionality, '_ncba.01')
TRB.v <- DNAStringSet(TRB.v)
TRB.v
writeXStringSet(TRB.v, './phylogen/fa/ncba.01/bovine_cattle.ncba01.TRB_V.fa')



#' ###  fetch TRG1 sequence of the new assembly (FWD in genome)

TRG1 <- readDNAStringSet('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/TRG/TRG_01_chr4_82.4_82.8Mb.fa')
TRG1.anno <- read_csv('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/TRG/TRG1_mapResults_clean.csv') %>% suppressMessages()

TRG1.v.anno <- TRG1.anno[grepl('TRGV', TRG1.anno$annotation),  ]
TRG1.v <- Views(TRG1[[1]], start = (TRG1.v.anno$start - 82400000 +1), width = TRG1.v.anno$length)

names(TRG1.v) <- paste0(TRG1.v.anno$annotation, '_', TRG1.v.anno$functionality, '_ncba.01')


TRG1.v<- DNAStringSet(TRG1.v)
TRG1.v
writeXStringSet(TRG1.v, './phylogen/fa/ncba.01/bovine_cattle.ncba01.TRG1_V.fa')



#' ###  fetch TRG2 sequence of the new assembly (rev in genome)

TRG2 <- readDNAStringSet('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/TRG/TRG_02_chr4_50.0_50.2Mb.fa')
TRG2.anno <- read_csv('./immuLocus/IG_TR/mapToGenome/geneDb/bt2/TRG/TRG2_mapResults_clean.csv')

TRG2.v.anno <- TRG2.anno[grepl('TRGV', TRG2.anno$annotation),  ]
TRG2.v <- Views(TRG2[[1]], start = (TRG2.v.anno$start - 50000000 +1), width = TRG2.v.anno$length) %>% reverseComplement()

names(TRG2.v) <- paste0(TRG2.v.anno$annotation, '_', TRG2.v.anno$functionality, '_ncba.01')

TRG2.v <- DNAStringSet(TRG2.v)
TRG2.v
writeXStringSet(TRG2.v, './phylogen/fa/ncba.01/bovine_cattle.ncba01.TRG2_V.fa')




