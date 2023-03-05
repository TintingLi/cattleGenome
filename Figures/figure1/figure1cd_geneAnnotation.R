


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(lessR))
#suppressPackageStartupMessages(library(karyoploteR))

#BiocManager::install('lessR')


#' - **figure 1c**  


layout(matrix(1))
par(mar=c(4,4,2,1))


TE <- data.frame(classes = c('totalRepeats', 'TE',  'LINE', 'SINE','LTR','SSR'),
				 ratio = c(47.30, 46.53, 28.70, 8.54,6.38, 0.08))
TE2 <- data.frame(classes = c('nonRepeats',  'LINE', 'SINE','LTR', 'other TE', 'tandem repeats'),
				  ratio = c(52.7,  28.70, 8.54,6.38,2.91, 0.41))

par(bg ='#CCFFFF')
par(bg = 'white')
#barplot(TE$ratio, ylim = c(0, 50), 
#		col = c('#003366', '#9999CC', '#66CCFF', '#CCCC33', '#FF6666', '#666666'))

PieChart(x=classes, y = ratio, data = TE2, lwd = 1,lty =1)
cols <- hcl.colors(length(unique(TE2$classes)))

PieChart(x=classes, y = ratio, data = TE2, fill = cols, 
		 values = 'off',init_angle = 45, clockwise = T)



#' - **figure 1d**  
#' 
gene.anno <- data.frame(classes = c('predicted', 'annotated',  'BUSCO_complete', 'expressed'),
						number = c(20288, 19794, 3943, 18193),
						ratio = c(1, 0.9757, 0.9608, 0.8967))

par(mar =c(4,4,2,1), oma = c(3,3,2,1))
barplot(gene.anno$ratio, col = c('#003366', '#9999CC', '#66CCFF', '#CCCC33'), horiz = T)



