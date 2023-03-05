
#' ---
#' title: "`mixcr` results visualization of pacBio RNA_Seq data"
#' author: "LiTingting(ting67@126.com)"
#' output:
#'   html_document:
#'     number_sections: false
#'     highlight: pygments
#'     theme: cerulean
#'     toc: yes
#' ---


suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(MetBrewer))

#options(width = 160)

setwd('~/share/litt/proj/cow/scripts/Figures/figure5/')
knitr::opts_knit$set(root.dir = "~/share/litt/proj/cow/scripts/Figures/figure5/")


#' - **figure 5b**  
mixcr.imgt<- read_tsv('./data/mixcr_summary_imgt.txt') %>% suppressMessages()
mixcr.imgt$sum <- apply(mixcr.imgt[, 2:5], 1, sum)
mixcr.imgt

mixcr.imgt.bp <- tibble(IGH= mixcr.imgt$sum[1:2], IGL = mixcr.imgt$sum[3:4],
						IGK= mixcr.imgt$sum[5:6], 
						TRA = c(mixcr.imgt$sum[9], as.integer(mixcr.imgt$sum[8] * (mixcr.imgt$sum[9] / mixcr.imgt$sum[7]))),
						TRD = c(mixcr.imgt$sum[10], as.integer(mixcr.imgt$sum[8] * (mixcr.imgt$sum[10] / mixcr.imgt$sum[7]))),
						TRB = mixcr.imgt$sum[11:12],
						TRG= mixcr.imgt$sum[13:14])
mixcr.imgt.bp
par(mar = c(2,2,1,0), oma = c(4,4,1,0))
barplot(mixcr.imgt.bp %>% as.matrix(), beside = T, 
		col = c('#62929a', '#859b6c'), ylim = c(0, 20000))

legend('topright', legend  = c('all', 'non-functional'), fill = c('#62929a', '#859b6c'))



#' - **figure 5c: CDRH3**  
#' 

cdr3 <- read_tsv('~/share/litt/proj/cow/results/pacbio_PBMC_RNA_seq/mixcr/imgt/clean/CDR3seqs_all.txt')
cdr3 <- cdr3[ !is.na(cdr3$nSeqCDR3) & 
			  	grepl('IGHV', cdr3$allVHitsWithScore) &
			  	grepl('IGHD', cdr3$allDHitsWithScore), ]
cdr3$CDR3Length  <- str_length(cdr3$nSeqCDR3)
cdr3
cdr3.ighd8_2 <- cdr3[grepl('IGHD8-2', cdr3$allDHitsWithScore) & !grepl(',', cdr3$allDHitsWithScore, ), ]
cdr3.non_ighd8_2 <- cdr3[!grepl('IGHD8-2', cdr3$allDHitsWithScore), ]
cdr3.ighv1_7 <- cdr3[grepl('IGHV1-7', cdr3$allVHitsWithScore) & !grepl(',', cdr3$allVHitsWithScore), ]
cdr3.nonighv1_7 <- cdr3[!grepl('IGHV1-7', cdr3$allVHitsWithScore), ]


palette(met.brewer('Tsimshian'))

par(oma = c(3,3,2,2))
plot(density(cdr3$CDR3Length/3,adjust = 2), xlim = c(0, 95), ylim = c(0, 0.08), col =1, lwd =3, axes =F)
lines(density(cdr3.ighd8_2$CDR3Length/3,adjust = 2), xlim = c(1, 95), col =2, lwd =3)
lines(density(cdr3.non_ighd8_2$CDR3Length/3,adjust = 2), xlim = c(1, 95), col =6, lwd =3)
lines(density(cdr3.nonighv1_7$CDR3Length/3,adjust = 2), xlim = c(1, 95), col =4, lwd =3)
lines(density(cdr3.ighv1_7$CDR3Length/3,adjust = 2), xlim = c(1, 95), col =7, lwd =3)

axis(side=1, at=seq(-5, 101, 10), line = -.5)
axis(side=2)

legend('topright', legend = c('All','IDHD8-2','IGHV1-7', 'Non IGHD8-2',  'Non IGHV1-7'), col =c(1,2,7, 6,4), lwd=4)




