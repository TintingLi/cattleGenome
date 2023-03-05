


suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RSkittleBrewer))
suppressPackageStartupMessages(library(RColorBrewer))
#suppressPackageStartupMessages(library(scales))
setwd('~/share/litt/proj/cow/results/')
knitr::opts_knit$set(root.dir = '~/share/litt/proj/cow/results/')




#' - **figure 5a: barplot of V gene numbers of three species**  
#' 
#' 
v_stat_total <- data.frame(hom=c(418, 219), mus=c(490, 329), bus=c(710, 351))
v_stat_IG <- data.frame(hom=c(279, 117), mus=c(334, 201), bus=c(201, 55))
v_stat_TR <- data.frame(hom=c(139, 102), mus=c(156, 128), bus=c(509, 296))

rownames(v_stat_total) <- c('total', 'F')
rownames(v_stat_IG) <- c('total', 'F')
rownames(v_stat_TR) <- c('total', 'F')

#+ fig.width=10, fig.height=5
layout(matrix(c(1,2,3), nrow = 1))

barplot(v_stat_total %>% as.matrix(), beside = T, col = c('#dadada', '#204969'), main = 'total')
barplot(v_stat_IG %>% as.matrix(), beside = T, col = c('#dadada', '#204969'), main = 'IG')
barplot(v_stat_TR %>% as.matrix(), beside = T, col = c('#dadada', '#204969'), main = 'TR')


layout(matrix(c(1,2), nrow = 1))
v.df <- data.frame(h_ig_v = c(117, 279-117), m_ig_v = c(201, 334-201), c_ig_v = c(55, 201-55),
				   h_tr_v =c(102, 139-102), m_tr_v = c(128, 156-128), c_tr_v = c(296, 509-296),
				   h_v = c(219, 418-219), m_v = c(329, 490 -329), c_v = c(351, 710-351))

tes <- barplot(v.df %>% as.matrix(), horiz = T, xlim = c(0, 800), col = c('#204969','#dadada'), las =2,
			   space = c(0.1, 0.1, 0.1, 0.618, 0.1, 0.1, 0.618, 0.1, 0.1))

axis(2, labels = F, at =c(-1, 10.936))

# 7.5 x 4.5
tes <- barplot(v.df %>% as.matrix(), horiz = F, ylim = c(0, 800), col = c('#204969','#dadada'), las =2,
			   space = c(0.1, 0.1, 0.1, 0.618, 0.1, 0.1, 0.618, 0.1, 0.1))

axis(1, labels = F, at =c(-1, 10.936), lwd.ticks = 0)

layout(matrix(1))






