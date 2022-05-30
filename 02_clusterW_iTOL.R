

#' **phylogenetic analysis  of V genes of three species**  
#' 
#'  + get all functional V gene sequences in `fa` format.  
#'  + multiple sequence alignment by `clusta`  
#'  + visualization of MSA results by `iTOL`  

#'
#' ### 1. merge all the functional V genes of cow assembly, human and mouse.  
#' 
#' `cat ncba.01/cattle.ncba01_V_Functionl_genes.fa human/human_V_functionalGenes.fa mouse/mouse_V_functionalGenes.fa > merged_V_functional.fa`
#' 


#' ### 2. `clusta` analysis  
#' <https://www.ebi.ac.uk/Tools/msa/clustalo/>  
#' up load the merged fa file and download the output.
#' 
#' ### 3. `iTOL` visualization.  
#'   + upload the clusta results into iTOL website.  
#'   + color annotations.
#'       - <https://github.com/BlakeRMills/MetBrewer>  
#'       - metbrewer$Austria for the color of 7 subgroups.  
#'       - metbrewer$Juarez for the color of 3 species.  




