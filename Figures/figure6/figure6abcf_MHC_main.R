#' ---
#' title: "gene structures of MHC"
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
suppressPackageStartupMessages(library(RSkittleBrewer))
suppressPackageStartupMessages(library(RColorBrewer))
options(width = 100)
setwd('~/share/litt/proj/cow/results/MHC/')
knitr::opts_knit$set(root.dir = '~/share/litt/proj/cow/results/MHC/')


annogenes <- read_csv('./geneAnno/mhc_bovin_gene_clean.csv') %>% suppressMessages()
annogenes <- filter(annogenes, identity > 60)
mhcgenes <- read.table('./boLA_align/19_markers/blast_Assembly/BoLA_gene_blast_cowAssembly_clean_merged.txt', 
					   col.names = c('genename','seqname', 'identity', 'alignmentLength',
					   			  'q.start','q.end','q.length', 's.start', 's.end', 's.length', 'evalue'))

bola.dra <- read.table('./boLA_align/BoLA_DRA/blast_Assembly/BoLA_DRA_gDNA_blast_cowAssembly.txt', 
					   col.names = c('genename','seqname', 'identity', 'alignmentLength',
					   			  'q.start','q.end','q.length', 's.start', 's.end', 's.length', 'evalue'))

mhcgenes <- rbind(mhcgenes, bola.dra)
mhcgenes <- mhcgenes[order(mhcgenes$s.start),]
write.table(mhcgenes, file = './boLA_align/final_merged_mhc_geneloci.txt', 
			quote = F, sep = '\t', row.names = F)

#' ### 1. figure_6a, upper panel,  plot the whole chromosome  
#' 

#+ fig.width=10, fig.height=5
par(oma =c(2,2,2,1))
plot(x=0,xlim=c(0, 54.2),ylim=c(0,3), type = 'n',yaxt = 'n', xaxt = 'n', 
	 main = 'MHC locus in the cow genome assembly',
	 xlab = "", ylab = '')

#add y label:
axis(1, at = c(0, 10, 20, 30, 40, 50, 54.2), cex=4, line = 1)

rate = 1000000
for (i in seq(1, dim(annogenes)[1])) {
#	print(i)
	rect(annogenes$start[i]/rate, 1, annogenes$end[i] / rate, 2, col=gray(0.8),lwd=0, border =gray(0.8) )
}


for (i in seq(1, dim(mhcgenes)[1])) {
#	print(i)
	rect(mhcgenes$s.start[i]/rate, 1, mhcgenes$s.end[i] / rate, 2, col='red',lwd=0.1, border = 'red')
	
}



#' ### 2. figure_6a, lower panel,  plot MHC Class I and II together.  
#' 
#+ fig.width=10, fig.height=3
par(mar =c(4,4,3,2))
plot(x=0,xlim=c(26,32),ylim=c(0,3), type = 'n',yaxt = 'n', xaxt = 'n', 
	 main = 'MHC locus in the cow genome assembly',
	 xlab = "", ylab = '')

#add y label:
axis(1, at = seq(26, 32, 0.5), cex=4)

rate = 1000000
for (i in seq(1, dim(annogenes)[1])) {
#	print(i)
	rect(annogenes$start[i]/rate, 1, annogenes$end[i] / rate, 2, col=gray(0.8),lwd=0, border =gray(0.8) )
}


for (i in seq(1, dim(mhcgenes)[1])) {
#	print(i)
	rect(mhcgenes$s.start[i]/rate, 1, mhcgenes$s.end[i] / rate, 2, col='red',lwd=0, border = 'red')
}



#' ### 3. figure_6b, plot MHC Class  I alone.  
#' 
#+ fig.width=10, fig.height=3

plot(x=0,xlim=c(28.7,30.5),ylim=c(0,3), type = 'n',yaxt = 'n', xaxt = 'n', 
	 main = 'MHC locus in the cow genome assembly',
	 xlab = "", ylab = '')

#add y label:
axis(1, at = seq(28.7,30.5, 0.1), cex=4)


#rect(0,Chr_id,Chr_length/ Rate,Chr_id + 0.6,col=gray(0.8),lwd=0)

rate = 1000000
for (i in seq(1, dim(annogenes)[1])) {
#	print(i)
	rect(annogenes$start[i]/rate, 1, annogenes$end[i] / rate, 2, col=gray(0.7),lwd=0, border =gray(0.8) )
}


for (i in seq(1, dim(mhcgenes)[1])) {
#	print(i)
	rect(mhcgenes$s.start[i]/rate, 1, mhcgenes$s.end[i] / rate, 2, col='red',lwd=0, border = 'red')
}


#' ### 4. figure_6c, plot MHC Class  II alone.  
#' 
#+ fig.width=8, fig.height=3

plot(x=0,xlim=c(25.7,27.5),ylim=c(0,3), type = 'n',yaxt = 'n', xaxt = 'n', 
	 main = 'MHC locus in the cow genome assembly',
	 xlab = "", ylab = '')

#add y label:
axis(1, at = seq(25.7,27.5, 0.1), cex=4)

rate = 1000000
for (i in seq(1, dim(annogenes)[1])) {
#	print(i)
	rect(annogenes$start[i]/rate, 1, annogenes$end[i] / rate, 2, col=gray(0.7),lwd=0, border =gray(0.8) )
}


for (i in seq(1, dim(mhcgenes)[1])) {
#	print(i)
	rect(mhcgenes$s.start[i]/rate, 1, mhcgenes$s.end[i] / rate, 2, col='red',lwd=0, border = 'red')
}



#' ### 5. figure_6f, plot gene locations on two haploids.  
#' 
#' haplotig locations

hp1 <- read.table('./haplotype/mergedhp1bam.tb.clean.txt', header = T) 
hp2 <- read.table('./haplotype/mergedhp2bam.tb.clean.txt', header = T)
head(hp1)

gene.hp1 <- read.table('./boLA_align/19_markers/blast_haplotype/geneMapToHP1_clean.txt', header = T) 
gene.hp2 <- read.table('./boLA_align/19_markers/blast_haplotype/geneMapToHP2_clean.txt', header = T)
head(gene.hp1)


#' **Add BoLA-DRA info**  
dra.hp1 <- read.table('./boLA_align/BoLA_DRA/blast_haplotype/BoLA_DRA_gDNA_mapto_hp1.txt')
colnames(dra.hp1) <- colnames(gene.hp1)
dra.hp1
dra.hp2 <- read.table('./boLA_align/BoLA_DRA/blast_haplotype/BoLA_DRA_gDNA_mapto_hp2.txt')
colnames(dra.hp2) <- colnames(gene.hp1)
dra.hp2

gene.hp1 <- rbind(gene.hp1, dra.hp1)
gene.hp2 <- rbind(gene.hp2, dra.hp2)

gene.hp1 <- left_join(gene.hp1,hp1, by = c('subject.acc.' = 'hp1ctg') )
gene.hp1$fstart <- 'NA'
gene.hp1$fstart[gene.hp1$strand == '+'] <- gene.hp1$s.start[gene.hp1$strand == '+'] +  gene.hp1$start[gene.hp1$strand == '+']
gene.hp1$fend[gene.hp1$strand == '+'] <- gene.hp1$s.end[gene.hp1$strand == '+'] + gene.hp1$start[gene.hp1$strand == '+']

gene.hp1$fend[gene.hp1$strand == '-'] <- gene.hp1$end[gene.hp1$strand == '-'] - gene.hp1$s.start[gene.hp1$strand == '-']
gene.hp1$fstart[gene.hp1$strand == '-'] <- gene.hp1$end[gene.hp1$strand == '-'] - gene.hp1$s.end[gene.hp1$strand == '-']
gene.hp1 <- gene.hp1[order(gene.hp1$fstart),]
head(gene.hp1)

min(c(1,2))


gene.hp2 <- left_join(gene.hp2,hp2, by = c('subject.acc.' = 'hp2ctg') )
gene.hp2$fstart <- 'NA'
gene.hp2$fstart[gene.hp2$strand == '+'] <- gene.hp2$s.start[gene.hp2$strand == '+'] + gene.hp2$start[gene.hp2$strand == '+']
gene.hp2$fend[gene.hp2$strand == '+'] <- gene.hp2$s.end[gene.hp2$strand == '+'] + gene.hp2$start[gene.hp2$strand == '+']

gene.hp2$fend[gene.hp2$strand == '-'] <- gene.hp2$end[gene.hp2$strand == '-'] - gene.hp2$s.start[gene.hp2$strand == '-']
gene.hp2$fstart[gene.hp2$strand == '-'] <- gene.hp2$end[gene.hp2$strand == '-'] - gene.hp2$s.end[gene.hp2$strand == '-']

gene.hp2 <- gene.hp2[order(gene.hp2$fstart),]
head(gene.hp2)

write_csv(gene.hp1, './haplotype/gene.hp1.csv')
write_csv(gene.hp2, './haplotype/gene.hp2.csv')


gene.hp1 <- read_csv('./haplotype/gene.hp1.csv') %>%suppressMessages()
gene.hp2 <- read_csv('./haplotype/gene.hp2.csv') %>%suppressMessages()


#+ fig.width =10, fig.height=5
par(mar =c(4,4,3,2))
plot(x=0,xlim=c(0,(32302131 -24918682)/rate),ylim=c(0,3), type = 'n',yaxt = 'n', xaxt = 'n', 
	 main = 'MHC locus in the cow genome assembly',
	 xlab = "", ylab = '')

#add y label:
axis(1, at = seq(0, 8, 0.5), cex=4)

for (i in seq(1, dim(hp2)[1])) {
	rect(hp2$start[i]/rate, 0.6, hp2$end[i]/rate, 1, col=gray(0.7),lwd=0, border =gray(0.8))
}

for (i in seq(1, dim(hp1)[1])) {
	rect(hp1$start[i]/rate, 2.4, hp1$end[i]/rate, 2, col=gray(0.7),lwd=0, border =gray(0.8))
}


for (i in seq(1, dim(gene.hp2)[1])) {
	tig <- gene.hp2$subject.acc.[i]
	tig.start <- hp2$start[hp2$hp2ctg == tig]
	rect((gene.hp2$s.start[i] + tig.start )/rate, 0.25, 
		 (gene.hp2$s.end[i] + tig.start )/rate, 0.5, col='black',lwd=1, border = 'black')
}


for (i in seq(1, dim(gene.hp1)[1])) {
	tig <- gene.hp1$subject.acc.[i]
	tig.start <- hp1$start[hp1$hp1ctg == tig]
	rect((gene.hp1$s.start[i] + tig.start )/rate, 2.75, 
		 (gene.hp1$s.end[i] + tig.start )/rate, 2.5, col= 'black',lwd=1, border = 'black')
}
