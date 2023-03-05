
#' **Figure 1B: genome coverage of ncba1.0 and gaps.**

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(GenomicRanges))
#suppressPackageStartupMessages(library(karyoploteR))
setwd('~/share/litt/proj/cow/results/')

knitr::opts_knit$set(root.dir = '~/share/litt/proj/cow/results/')



#' ## genome length:
new.genomeLen <- read_csv('./metadata/NewAssembly//assembly_chrLength.csv') %>%
				 filter(!grepl('Scaffold', seqname)) %>% suppressMessages()
new.genomeLen$accum.len <- c(1, sapply(1:29, function(x){sum(new.genomeLen$length[1:x])}))
new.genomeLen

bos9.genomeLen <- read_delim('./metadata/ARS.UCD1.2/bosTau9.chrom.sizes.clean.txt', delim = '\t',
							 col_names  = c( 'seqname','length'))%>% suppressMessages()
bos9.genomeLen$accum.len <- c(1, sapply(1:29, function(x){sum(bos9.genomeLen$length[1:x])}))
bos9.genomeLen

#' ## gap loci of the new assembly and bostau9

assembly.gaploci <- read_csv('./metadata/NewAssembly//cowNewAssembly_N_gap_locations.csv') %>% 
					select(any_of(c('seqnames', 'start'))) %>% filter(!grepl('Scaffold', seqnames)) %>% suppressMessages()
assembly.gaploci
bos9.gaploci <- read_csv('./metadata/ARS.UCD1.2/ARS_UCD1.2_N_gap_locations.csv') %>%
					select(any_of(c('seqnames', 'start')))  %>% filter(!grepl('chrUn_', seqnames)) %>% suppressMessages()
bos9.gaploci

#' adjust the bostau9 loci according to the new assembly.
bos9.gaploci$adj.start <- sapply(1:length(bos9.gaploci$start), function(x){
	chrname = bos9.gaploci$seqnames[x]
	return(bos9.gaploci$start[x] * 
		   (new.genomeLen$length[new.genomeLen$seqname ==chrname] / 
		    	bos9.genomeLen$length[bos9.genomeLen$seqname == chrname])  )
})


assembly.gaploci$accum.loci <- sapply(1:dim(assembly.gaploci)[1], function(x){
	return(assembly.gaploci$start[x] + 
		   	new.genomeLen$accum.len[new.genomeLen$seqname == assembly.gaploci$seqnames[x]])
})


bos9.gaploci$accum.loci <- sapply(1:dim(bos9.gaploci)[1], function(x){
	return(bos9.gaploci$start[x] + 
		   	bos9.genomeLen$accum.len[bos9.genomeLen$seqname == bos9.gaploci$seqnames[x]])
})


bos9.gaploci$accum.adj.loci <- sapply(1:dim(bos9.gaploci)[1], function(x){
	return(bos9.gaploci$adj.start[x] + 
		   	new.genomeLen$accum.len[new.genomeLen$seqname == bos9.gaploci$seqnames[x]])
})

#' ## genome coverage info


assembly_cov_50k <- read_csv('./metadata/depthCovGCcontent/mosdepth/assembly/cow_assembly_coverage_50kb_q10.csv') %>% suppressMessages()
assembly_cov_50k

#' ## visualiztion.
#layout(matrix(1))
#+ fig.width=18, fig.height=3, fig.align = 'center'
layout(matrix(c(1,1,2,2,3,3,4,4,4), nrow =9, byrow = T))
par(mar =c(1,4,1,2), oma = c(1,2,1,1))

plot(x=0,xlim=c(0,dim(assembly_cov_50k)[1]),ylim=c(0,300), type = 'n', axes =F, xlab = "", ylab = '')
rect(1:dim(assembly_cov_50k)[1], 0.1, (1:dim(assembly_cov_50k)[1]) +1, assembly_cov_50k$ont_50kb_cov, lwd=1, col = '#FFCCCC', border = '#FFCCCC')
axis(2, at = c(0,300), line = -3)

#plot.new()
plot(x=0,xlim=c(0,dim(assembly_cov_50k)[1]),ylim=c(0,50), type = 'n', axes =F, xlab = "", ylab = '')
rect(1:dim(assembly_cov_50k)[1], 0.1, (1:dim(assembly_cov_50k)[1]) +1, assembly_cov_50k$ccs_50kb_cov, lwd=1, col = '#00CCCC', border = '#00CCCC')
axis(2, at = c(0,50), line = -3)

plot(x=0,xlim=c(0,dim(assembly_cov_50k)[1]),ylim=c(0,200), type = 'n', axes =F, xlab = "", ylab = '')
rect(1:dim(assembly_cov_50k)[1], 0.1, (1:dim(assembly_cov_50k)[1]) +1, assembly_cov_50k$ngs_50kb_csv, lwd=1, col = '#009966', border = '#009966')
axis(2, at = c(0,200), line = -3)


plot(x=0,xlim=c(0,sum(new.genomeLen$length)),ylim=c(0, 0.8), type = 'n', axes =F, 
	 xlab = "", ylab = '')
axis(1, at = c(new.genomeLen$accum.len, sum(new.genomeLen$length)),labels = F,  cex=2,2)

rect(bos9.gaploci$accum.adj.loci, 0.1, bos9.gaploci$accum.adj.loci +1, 0.4, lwd=1)
rect(assembly.gaploci$accum.loci, 0.5, assembly.gaploci$accum.loci +1, 0.8, lwd=1)

layout(matrix(1))



