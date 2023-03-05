
#' **Figure 1A: Global picture of NCBA1.0 assembly.**


################################################################################
#
# **Plot the Global picture of the new assembly** 
#
#   1. contigs loci  
#   2. segmental duplicates  
#   3. filling gaps
#   4. placed scaffolds
#   5. satDNAs
#   5. telomeres
#
################################################################################


suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RSkittleBrewer))
suppressPackageStartupMessages(library(RColorBrewer))
#suppressPackageStartupMessages(library(scales))
setwd('~/share/litt/proj/cow/results/')
knitr::opts_knit$set(root.dir = '~/share/litt/proj/cow/results/')

seqlength <- read_tsv('./metadata/NewAssembly/chr_length_gcContent.txt',
					  col_names = c('seqname', 'length', 'gcContent'))

# ctglocus ----
ctglocus <- read_tsv('./contigMaptoGenome/ctgMapDetails_update.txt')
seqnames <- ctglocus$seqname %>% unique() %>% sort()
ctglocus$width <- ctglocus$seqEnd - ctglocus$seqStart
ctglocus %<>% arrange(width)
ctglocus$color <- 1:184

# segmental duplicates ----
segdups <- read_tsv('./segdups/sedef/newAssemblySoftSegDup_sedeg.clean',
					col_names = c('seq1', 'start1', 'end1',
								  'seq2', 'start2', 'end2',
								  'sdname','score','strand1','strand2',
								  'max_len','aln_len','comment',
								  'indel_a','indel_b', 'alnB',
								  'matchB','mismatchB','transitionsB',
								  'transversions', 'fracMatch',
								  'fracMatchIndel', 'jck', 'K2K',
								  'aln_gaps', 'uppercaseA', 'uppercaseB',
								  'uppercaseMatches','aln_matches', 'aln_mismatches',
								  'aln_gaps2', 'aln_gap_bases', 'cigar', 'filter_score')) %>% suppressMessages()

segdups <- filter(segdups, aln_len > 5000)

segcords <- rbind(segdups[,1:3] %>% set_colnames(c('seqname', 'seqStart', 'seqEnd')),
				  segdups[,4:6]%>% set_colnames(c('seqname', 'seqStart', 'seqEnd'))) 
segcords <- filter(segcords, !grepl('Scaffold', x = segcords$seqname))

# gaploci ---------
gaploci.final <- read_tsv('./FillingGaps/ars_ucd1.2_chr_ngaps/ARS_UCD1.2_Ngap_loci_in_Assebmly_final.txt',
						  col_names = c('seqname', 'start', 'end', 'Gap_ars_ucd1.2', 
						  			  'ctg_seqname', 'ctg_start','ctg_end','ctg')) %>% suppressMessages()
# scaffold loci ------
scaffold.loci.final <- 
	read_tsv('./FillingGaps/Ars.ucd1.2.scaffolds/ARS_UCD1.2_scafold_loci_in_Assebmly_final.txt',
			 col_names = c('seqnames', 'start', 'end', 'strand', 'qwidth', 'width', 'njuc', 
			 			  'scaffold', 'scaffoldLength', 'locus_chr', 'locus_chr_ctg_start', 
			 			  'locus_chr_ctg_end', 'locus_ctg')) %>% suppressMessages()
scaffold.loci.final

# satDNA ------

satDNA.loci <- read_tsv('./centromere/bt2/satLociMerge.tsv', col_names = c('seqname', 'start','end'))

satDNA.loci <- filter(satDNA.loci, start != 'na')
satDNA.loci$start <- as.integer(satDNA.loci$start)
satDNA.loci$end <- as.integer(satDNA.loci$end)




#' **plots**
#' 
par(mar =c(4,4,3,2))
#+ fig.width=10, fig.height=10
Rate <- 1e6
plot(x=0,xlim=c(0,170),ylim=c(0,60), type = 'n',yaxt = 'n', xaxt = 'n', 
	 main = 'contig locus in the cow genome assembly',
	 xlab = "", ylab = '')

#add y label:
axis(1, at = c(0, 50, 100, 150, 180), 
	 labels = c('0','50M','100M','150M','180M'), cex=2)
axis(2, at = seq(0.5, 60, by =2), labels = seqnames[1:30], 
	 tick = F, las =2, line = -0.7, cex=2)

# **1. chrs and contigs **  
Chr_id <- 0
for(i in seqnames[1:30]){
	
	Chr_length <- seqlength$length[seqlength$seqname == i] 
#	print(i)
#	print(Chr_length)
	
	#plot the whole chr:
	rect(0,Chr_id,Chr_length/ Rate,Chr_id + 1,col=gray(0.4),lwd=0)
	
	selectcontig <- filter(ctglocus, seqname == i)
	
	if(dim(selectcontig)[1] == 1){
		#one chromosome contains only one contig:
		X1 <- selectcontig[1,'seqStart'] / Rate
		X2 <- selectcontig[1, 'seqEnd'] / Rate
		Y1 <- Chr_id 
		Y2 <- Chr_id + 1
		rect(X1,Y1,X2,Y2,col=brewer.pal(9,"YlGnBu")[9],lwd=1)
	}else {
		for(i in seq(1,dim(selectcontig)[1])) {
			X1 <- selectcontig[i,'seqStart'] / Rate
			X2 <- selectcontig[i, 'seqEnd'] / Rate
			Y1 <- Chr_id 
			Y2 <- Chr_id + 1
			tmp <- selectcontig[i,]
			
			if(tmp$width > 50000000){
				scol <- brewer.pal(9,"YlGnBu")[7] 
			}else if ( (tmp$width < 50000000) & (tmp$width > 20000000)){
				scol <- brewer.pal(9,"YlGnBu")[6] 
			}else if ((tmp$width < 20000000) & (tmp$width > 2000000)){
				scol <- brewer.pal(9,"YlGnBu")[4] 
			}else if ((tmp$width > 2000000) & (tmp$width > 500000)){
				scol <- brewer.pal(9,"YlGnBu")[2] 
			}else {
				scol <- brewer.pal(9,"YlGnBu")[1] 
			}
			
			rect(X1,Y1,X2,Y2,col=scol,lwd=1)
		}}
	Chr_id <- Chr_id + 2
}



# **2. segmental duplicates **  

Chr_id <- 0
for(i in seqnames[1:30]){
	
	tmp <- filter(segcords, seqname == i)
	for(i in seq(1,dim(tmp)[1])) {
		X1 <- as.numeric(tmp[i,'seqStart']) / Rate
		X2 <- as.numeric(tmp[i, 'seqEnd'] )/ Rate
		Y1 <- Chr_id +1
		Y2 <- Chr_id + 1.3
		rect(X1,Y1,X2,Y2,lwd=1)
	}
	Chr_id <- Chr_id + 2
}



# **3. filled gaps**.  

Chr_id <- 0
for(i in seqnames[1:30]){
	
	tmp <- filter(gaploci.final, seqname == i)
	
	for(i in seq(1,dim(tmp)[1])) {
		X1 <- tmp[i,'start'] / Rate
		Y1 <- Chr_id
		points(X1, Y1, pch=17, col = '#FF0000', cex = 0.6)
	}
	Chr_id <- Chr_id + 2
}

# **4. placed scaffolds**.  

Chr_id <- 0
for(i in seqnames[1:30]){
	
	tmp <- filter(scaffold.loci.final, seqnames == i)
	
	for(i in seq(1,dim(tmp)[1])) {
		X1 <- c(tmp[i,'start'] / Rate, tmp[i,'end'] / Rate)
		Y1 <- Chr_id
		#    points(X1, c(Y1+0.8, Y1+0.8),  col =  "#FFDC14", cex = 0.6, type = 'l', lwd =4)
		rect(tmp[i,'start'] / Rate, Chr_id+ 0.5, tmp[i,'end'] / Rate, Chr_id+ 1,col =  "#FFDC14", lwd =0 )
	}
	Chr_id <- Chr_id + 2
}

# **5. satDNA**
#
Chr_id <- 0
for(i in seqnames[1:30]){
	
	tmp <- filter(satDNA.loci, seqname == i)
	if(dim(tmp)[1] != 1){
		Chr_id <- Chr_id + 2
		next()
	} else{
		
		rect(tmp$start / Rate, Chr_id + 0.25, tmp$end / Rate, Chr_id + 0.75, col = '#FF6600', lwd=0)
		Chr_id <- Chr_id + 2
	}
}


# **6. telomeres**. 
# 

teloChrs <- c('chr06', 'chr17', 'chr18', 'chr20', 'chr22', 'chr26', 'chr28')

Chr_id <- 0
for(i in seqnames[1:30]){
	
	
	if( !(i %in% teloChrs) ){
		Chr_id <- Chr_id + 2
		next()
	} else{
		tmp <- filter(seqlength, seqname == i)
		rect(tmp$length / Rate, Chr_id, tmp$length / Rate + 0.5, Chr_id + 1, col = '#9933FF', lwd=0)
		Chr_id <- Chr_id + 2
	}
}










