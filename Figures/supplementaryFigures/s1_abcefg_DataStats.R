#' ---
#' title: "Statistics of sequenced data of the Cow Genome Project"
#' author: "LiTingting(ting67@126.com)"
#' output: 
#'   html_document:
#'     number_sections: false
#'     highlight: pygments
#'     theme: cerulean
#'     toc: yes
#' ---


suppressPackageStartupMessages(library(tidyverse) )
library(magrittr) %>% suppressPackageStartupMessages()
library(RSkittleBrewer)%>% suppressPackageStartupMessages()
library(readxl)%>% suppressPackageStartupMessages()
library(lubridate)%>% suppressPackageStartupMessages()
library(plotrix)%>% suppressPackageStartupMessages()
setwd('~/share/litt/proj/cow/')
knitr::opts_knit$set(root.dir = '~/share/litt/proj/cow/')




#' ### 1. figure_s1a, readLength distribution per cell.   
#'   
#'  - ONT ultra-long reads  
#'     - batch1: 15 cells, total 237.8Gb, N50: 70.4Kb  
#'     - batch2: 15 cells, total 261.1Gb,   
#'     - total: 30 cells, total 499.0Gb.  
#' - pacbio hifi reads    
#'     - total 17 cells.  
#' 
ONT_ultra_long <- read_excel("./results/metadata/stats/ONT_ultra_long_merged.xlsx")
ONT_ultra_long$library_ID
ONT_ultra_long$TotalPassReadsBases[1:15] %>% sum()
ONT_ultra_long$TotalPassReadsNumbers[1:15] %>% sum()
ONT_ultra_long$TotalPassReadsBases[16:30] %>% sum()
ONT_ultra_long$TotalPassReadsNumbers[16:30] %>% sum()
sum(ONT_ultra_long$TotalPassReadsBases)
mean(ONT_ultra_long$passReadsN50Length)

ccs <- read_excel('./results/metadata/stats/ccs_stat.xlsx')
ccs$Movies

#' N50 length :
palette(c(RSkittleBrewer('wildberry'),RSkittleBrewer('tropical') ))
midpoints <- 
	barplot(c(ONT_ultra_long$passReadsN50Length, ccs$SubReadsN50Length),
			ylim = c(0, 130000), las=2,
			col = rep(c(3,2), times = c(30,17)))

# N50 of CCS  
abline(h = 11757, lty = 'dashed', lwd =1.5, col = 'gray')
# N50 of ONT reads
lines(x = c(0,37), y = c(70386, 70386), lty = 'dashed', lwd =1.5, col = 'gray')

#' ### 2. figure_s1b, read yield per cell.   
#' total read yields:

layout(matrix(c(1)))
palette(c(RSkittleBrewer('wildberry'),RSkittleBrewer('tropical') ))
barplot(c(ONT_ultra_long$TotalPassReadsBases, ccs$totalSubReadBases),
		ylim = c(0, 40000000000),las=2,
		col = rep(c(3,2), times = c(30,17)))


#' ### 3. figure_s1c, read length density plot.   

#ontReads <- read.table('./results/stats/readLength/ont.bam.readLength.stat')
ontReads <- read.table('./results/metadata/stats/readLength/tmp/ont_merged.readlength.stat')

names(ontReads) <- c('readNum', 'readLength')
ontReads <- as_tibble(ontReads)
ontReads
readLengthDis <- tibble(length = seq(1000, 900000, 1000))

readLengthDis$readNum <- sapply(1:length(readLengthDis$length), function(x){
	if(x ==1){ return(
		sum(ontReads$readNum[ontReads$readLength <= readLengthDis$length[1]]))
	}else{return(
		sum(ontReads$readNum[(ontReads$readLength <= readLengthDis$length[x]) &
							 	(ontReads$readLength > readLengthDis$length[x-1])])
	)}
})

#+ fig.width=5, fig.height=5, fig.align = 'center'
#plot(readLengthDis$length[1:200], readLengthDis$readNum[1:200], type = 'h')
#axis.break(2, c(3e5))
par(oma = c(4,4,2,1))
#plot(readLengthDis$length[1:200], readLengthDis$readNum[1:200], type = 'h',  col ='#352B73', lwd =2, las=2)
plot(readLengthDis$length, readLengthDis$readNum, type = 'h',  col ='#352B73', lwd =2, las=2)

#axis.break(2, c(1.5e5))


#' ### 4. figure s1e, contig length accumulation 
ctgLength <- read_tsv('./results/metadata/NewAssembly/contigLengthBED.txt') %>% arrange(-end)
ctgLength$acuum <- sapply(1:length(ctgLength$end), function(x){return(sum(ctgLength$end[1:x]))})
#+ fig.width=5, fig.height=5, fig.align = 'center'
layout(matrix(1))
plot(ctgLength$acuum, ctgLength$end, type = 'l', lwd =3, col = 6)



#' ### 5. figure s1f, global coverage of the genome.
#' 
layout(matrix(1))
mosresultsfiles <- list.files('./results/metadata/depthCovGCcontent/mosdepth/assembly/mosdepth_q10', 
								pattern = ".*50k.*global.dist.txt",	 full.names = T)


mosresults <- lapply(mosresultsfiles, function(x) {read.table(x) %>% 
		as_tibble %>% filter(V1=='total')})

names(mosresults) <- basename(mosresultsfiles) %>% 
	gsub('.mosdepth.global.dist.txt', '', x = .) 

palette(RSkittleBrewer("M&M"))
#+  fig.width = 6, fig.height=6, fig.align= 'center'  
#par(mar =c(1,1,1,1), oma =c(2,2,2,2))
plot(mosresults[[1]]$V2, mosresults[[1]]$V3, xlim = c(0, 250), ylim=c(0,1),
	 type = "l", lwd=4, col='#FF6666', las=1, ann =F,
	 xlab = 'genome coverage',
	 ylab = 'proportion of genome at coverage'#, main = 'Coverage analysis'
)
lines(mosresults[[2]]$V2, mosresults[[2]]$V3, 
	  xlim = c(0, 10), type = "l", lwd=4, col = '#663333')
lines(mosresults[[3]]$V2, mosresults[[3]]$V3, 
	  xlim = c(0, 10), type = "l", lwd=4, col ='#003399')

legend(x = 200, y=1, 
	   legend = c('ccs', 'ont', 'ngs') , 
	   lwd=2, col = c('#FF6666',  '#003399','#663333'), cex=1, bty = 'n' )


#' only keep CCS and ONT:
#+  fig.width = 6, fig.height=6, fig.align= 'center'  
#par(mar =c(1,1,1,1), oma =c(2,2,2,2))
plot(mosresults[[1]]$V2, mosresults[[1]]$V3, xlim = c(0, 250), ylim=c(0,1),
	 type = "l", lwd=4, col='#FF6666', las=1, ann =F,
	 xlab = 'genome coverage',
	 ylab = 'proportion of genome at coverage'#, main = 'Coverage analysis'
)
#lines(mosresults[[2]]$V2, mosresults[[2]]$V3, 
#	  xlim = c(0, 10), type = "l", lwd=4, col = '#663333')
lines(mosresults[[3]]$V2, mosresults[[3]]$V3, 
	  xlim = c(0, 10), type = "l", lwd=4, col ='#003399')

legend(x = 160, y=1, 
	   legend = c('ccs', 'ont') , 
	   lwd=2, col = c('#FF6666',  '#003399'), cex=1, bty = 'n' )



#' ### 6 figure s1g, correlations between gcContent and genome coverage.
#' 

ctg500bpClean <- read_csv('./results/metadata/depthCovGCcontent/ctg500bpclean.csv')
ctg500bpClean

mean_depth <- ctg500bpClean %>% summarise(ccs_mean_depth = sum(n*ccs_mean) /sum(n),
										  ont_mean_depth = sum(n*ont_mean) /sum(n),
										  ngs_mean_depth = sum(n*ngs_mean) /sum(n))

#+  fig.width = 6, fig.height=6, fig.align= 'center'  
#palette(RSkittleBrewer("tropical"))
palette(RSkittleBrewer("M&M"))
par(mar = c(4,4,4,4), oma = c(2,2,2,2))
plot(ctg500bpClean$gr, ctg500bpClean$GCratio ,
	 ylim=c(0, 0.1), xlim =c(0,100),xaxt='n', type ='h', lwd=2,
	 col = '#66CC33', las=1,
	 xlab = '%GC of 500bp bins',
	 ylab = 'fraction of bins of %GC',
	 main ='GC contents VS mapping coverage')
axis(1, at = seq(0,100, 10))

par(new =T)
plot(ctg500bpClean$gr, log10(ctg500bpClean$ont_mean/mean_depth$ont_mean_depth) , 
	 col ='#003399',  lwd=4, ylim = c(-3,3), xlim = c(0,100),
	 type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")

lines(ctg500bpClean$gr, log10(ctg500bpClean$ngs_mean/mean_depth$ngs_mean_depth) ,
	  lwd=4, col = '#663333')
lines(ctg500bpClean$gr, log10(ctg500bpClean$ccs_mean/mean_depth$ccs_mean_depth) ,
	  lwd=4, col ='#FF6666')

axis(4, at = seq(-3,3,1))
mtext('log10(coverage/(mean depth))', side = 4, line=3)

legend(x = 50, y=3, 
	   legend = c('ccs', 'ont', 'ngs'), 
	   lwd=2,  col = c('#FF6666',  '#003399','#663333'), cex=1, bty = 'n' )

