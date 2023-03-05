#' ---
#' title: "haplotypes of MHC"
#' author: "LiTingting(ting67@126.com)"
#' output: 
#'   html_document:
#'     number_sections: false
#'     highlight: pygments
#'     theme: cerulean
#'     toc: yes
#' ---

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Biostrings))
options(width = 100)
setwd('~/share/litt/proj/cow/results')
knitr::opts_knit$set(root.dir = '~/share/litt/proj/cow/results')




#' # 1. figure_6e,  **contig Length distribution**  

mhc <- readDNAStringSet('../data/nextomics/cow_anno_final/MHC_Haplotype/result/chr23.mhc.fasta')
hp1 <- readDNAStringSet('../data/nextomics/cow_anno_final/MHC_Haplotype/result/niu.HP1.contigs.fasta')
hp2 <- readDNAStringSet('../data/nextomics/cow_anno_final/MHC_Haplotype/result/niu.HP2.contigs.fasta')


layout(matrix(c(1,2)))
par(mar = c(4,2,2,1), mar= c(4,2,2,1))
barplot(width(hp1) %>% sort(decreasing = T), col = '#106898', las=1)
barplot(width(hp2)  %>% sort(decreasing = T), col = '#003460', las=1)



#' # 2. **aligned hp contigs back to the MHC assembly**    
#' hp1: 
hp1ctgLength <- read_tsv('./MHC/haplotype/hp1_ctgLength.txt', col_names = c('hp1ctg', 'hp1ctgLength')) %>% suppressMessages()
hp1bam <- readGAlignments('./MHC/haplotype/hp1_mapToMHC.bam', use.names = T)

tes <- reduce(granges(hp1bam))
tes$ctgs <- "tig1"

mergedhp1bam <- GRanges()
for (i in unique(names(hp1bam))) {
  tmp <- hp1bam[names(hp1bam) == i] %>% GRanges() %>% reduce(min.gapwidth=1000)
  tmp$hp1ctg <- i
  mergedhp1bam <- append(mergedhp1bam, tmp)
  
}

mergedhp1bam.tb <- as.data.frame(mergedhp1bam) %>% as_tibble()
mergedhp1bam.tb <- left_join(mergedhp1bam.tb, hp1ctgLength)
write_tsv(mergedhp1bam.tb, './MHC/haplotype/mergedhp1bam.tb.txt')

#' hp2: 
hp2ctgLength <- read_tsv('./MHC/haplotype/hp2_ctgLength.txt', col_names = c('hp2ctg', 'hp2ctgLength')) %>% suppressMessages()
hp2bam <- readGAlignments('./MHC/haplotype/hp2_mapToMHC.bam', use.names = T)

tes <- reduce(granges(hp2bam))
tes$ctgs <- "tig1"

mergedhp2bam <- GRanges()
for (i in unique(names(hp2bam))) {
  tmp <- hp2bam[names(hp2bam) == i] %>% GRanges() %>% reduce(min.gapwidth=1000)
  tmp$hp2ctg <- i
  mergedhp2bam <- append(mergedhp2bam, tmp)
  
}

mergedhp2bam.tb <- as.data.frame(mergedhp2bam) %>% as_tibble()
mergedhp2bam.tb <- left_join(mergedhp2bam.tb, hp2ctgLength)

write_tsv(mergedhp2bam.tb, './MHC/haplotype/mergedhp2bam.tb.txt')



#' # 3. figure_6d    
#' 

hp1loci <- read_tsv('./MHC/haplotype/mergedhp1bam.tb.clean.txt')
hp1loci$ystart <- sapply(1:dim(hp1loci)[1], function(x){return(sum(hp1loci$width[0 : (x-1)]))})
hp1loci$yend <- sapply(1:dim(hp1loci)[1], function(x){return(sum(hp1loci$width[1:x]))})

hp2loci <- read_tsv('./MHC/haplotype/mergedhp2bam.tb.clean.txt')
hp2loci$ystart <- sapply(1:dim(hp2loci)[1], function(x){return(sum(hp2loci$width[0 : (x-1)]))})
hp2loci$yend <- sapply(1:dim(hp2loci)[1], function(x){return(sum(hp2loci$width[1:x]))})

tail(hp1loci)
mhcLength <- 32302131 - 24918682
mhcLength #7383449

layout(matrix(1))

#+ fig.width=5, fig.height =5
plot.new()
par(mar =c(4,4,3,2))
plot(x=0,xlim=c(0,7383449),ylim=c(0,7100000), type = 'n',yaxt = 'n', xaxt = 'n', 
#     main = 'contig locus in the cow genome assembly',
     xlab = "", ylab = '')

#add y label:
axis(1, at = c(0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000), 
     labels = c('0','1M','2M','3M','4M','5M','6M','7M'), cex=2)

for (i in 1:dim(hp1loci)) {
  tmp <- hp1loci[i, ]
  abline(h = tmp$ystart, lty = 'dotted', lwd=0.5,  col ='grey')
#  abline(v = tmp$yend)
  if(tmp$strand == '+'){
    lines(x = c(hp1loci$start[i], hp1loci$end[i]),
        y = c(hp1loci$ystart[i], hp1loci$yend[i]), type = 'l' , col = '#5AA4AE', lwd= 3)
    points(hp1loci$start[i], hp1loci$ystart[i], pch= 19, col = '#5AA4AE')
    points(hp1loci$end[i], hp1loci$yend[i], pch= 19, col = '#5AA4AE')
  }else{
    lines(x = c(hp1loci$start[i], hp1loci$end[i]),
          y = c( hp1loci$yend[i], hp1loci$ystart[i]), type = 'l' , col = '#5AA4AE', lwd= 3)
    points(hp1loci$start[i], hp1loci$yend[i], pch= 19, col = '#5AA4AE')
    points(hp1loci$end[i], hp1loci$ystart[i], pch= 19, col = '#5AA4AE') 
    
  }
}

for (i in 1:dim(hp2loci)) {
  tmp <- hp2loci[i, ]
  abline(h = tmp$ystart, lty = 'dotted', lwd=0.5,  col ='grey')
  if(tmp$strand == '+'){
    lines(x = c(hp2loci$start[i], hp2loci$end[i]),
          y = c(hp2loci$ystart[i], hp2loci$yend[i]), type = 'l' , col = '#FEC85E', lwd= 3)
    points(hp2loci$start[i], hp2loci$ystart[i], pch= 19, col = '#FEC85E')
    points(hp2loci$end[i], hp2loci$yend[i], pch= 19, col = '#FEC85E')
  }else{
    lines(x = c(hp2loci$start[i], hp2loci$end[i]),
          y = c( hp2loci$yend[i], hp2loci$ystart[i]), type = 'l' , col = '#FEC85E', lwd= 3)
    points(hp2loci$start[i], hp2loci$yend[i], pch= 19, col = '#FEC85E')
    points(hp2loci$end[i], hp2loci$ystart[i], pch= 19, col = '#FEC85E') 
    
  }
}



for (i in 1:dim(hp2loci)) {
  tmp <- hp2loci[i, ]
  abline(h = tmp$ystart, lty = 'dotted', lwd=0.5,  col ='grey')
  if(tmp$strand == '+'){
    lines(x = c(hp2loci$start[i], hp2loci$end[i]),
          y = c(hp2loci$ystart[i], hp2loci$yend[i]), type = 'l' , col = '#887657', lwd= 3)
    points(hp2loci$start[i], hp2loci$ystart[i], pch= 19, col = '#887657')
    points(hp2loci$end[i], hp2loci$yend[i], pch= 19, col = '#887657')
  }else{
    lines(x = c(hp2loci$start[i], hp2loci$end[i]),
          y = c( hp2loci$yend[i], hp2loci$ystart[i]), type = 'l' , col = '#7A675E', lwd= 3)
    points(hp2loci$start[i], hp2loci$yend[i], pch= 19, col = '#7A675E')
    points(hp2loci$end[i], hp2loci$ystart[i], pch= 19, col = '#7A675E') 
    
  }
}


for (i in c(0,1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000 )) {
  abline(v=i, lty ='dotted', lwd =0.5, col ='grey')
}














