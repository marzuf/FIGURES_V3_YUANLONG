
options(scipen=100)

# Rscript cmp_meanCorr_permut_v24_ratioDown.R

script_name <- "cmp_meanCorr_permut_v24_ratioDown.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildData <- TRUE

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

outFolder <- file.path("CMP_MEANCORR_PERMUT_V24_RATIODOWN")
dir.create(outFolder, recursive = TRUE)

runFolder <- file.path("..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA")

permutFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")

plotMargin <- c(1,2,1,1)
legTitle <- "FCC ranges:"
fractBarSubTitle <- "AUC ratios:\n"
fractBarTitle_main <- "Fold-change concordance scores - PERMUT DATA"

ggsci_pal <- "lancet"
ggsci_subpal <- ""

cumsumCurve_col <- "darkgrey"


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "svg"

myHeight <- 7
myWidth <- 10
myHeightGG <- 7
myWidthGG <- 10

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")


all_hicds <- list.files(file.path(permutFolder))
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(permutFolder, x)))


if(buildData) {
  all_permut_fcc <- foreach(hicds = all_hicds) %do%{
    all_exprds_fcc <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
     
      all_data_v2 <- get(load(file.path("RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V2", hicds, exprds, "ds_all_permut.Rdata")))
      all_ratioDown_v2a <- unlist(lapply(all_data_v2, function(x) x[["ratioDown_right"]]))
      all_ratioDown_v2b <- unlist(lapply(all_data_v2, function(x) x[["ratioDown_left"]]))
      all_ratioDown_v2 <- c(all_ratioDown_v2a, all_ratioDown_v2b)
      
#      all_data_v3 <- get(load(file.path("RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V3", hicds, exprds, "ds_all_permut.Rdata")))
#      all_ratioDown_v3a <- unlist(lapply(all_data_v3, function(x) x[["ratioDown_right"]]))
#      all_ratioDown_v3b <- unlist(lapply(all_data_v3, function(x) x[["ratioDown_left"]]))
#      all_ratioDown_v3 <- c(all_ratioDown_v3a, all_ratioDown_v3b)
      
      all_data_v2_either <- get(load(file.path("RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V2_EITHER", hicds, exprds, "ds_all_permut.Rdata")))
      all_ratioDown_v2_either <- unlist(lapply(all_data_v2_either, function(x) x[["ratioDown"]]))
      
#      all_data_v3_either <- get(load(file.path("RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V3_EITHER", hicds, exprds, "ds_all_permut.Rdata")))
#      all_ratioDown_v3_either <- unlist(lapply(all_data_v3_either, function(x) x[["ratioDown"]]))
      
      all_data_v4_either <- get(load(file.path("RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V4_EITHER", hicds, exprds, "ds_all_permut.Rdata")))
      all_ratioDown_v4_either <- unlist(lapply(all_data_v4_either, function(x) x[["ratioDown"]]))
      
      
      list(
        ratioDown_meanRL = all_ratioDown_v2,
#        ratioDown_meanRL_withTAD = all_ratioDown_v3,
        ratioDown_RorL = all_ratioDown_v2_either,
#        ratioDown_RorL_withTAD = all_ratioDown_v3_either,
        ratioDown_RorL_replace = all_ratioDown_v4_either
      )
      
    }
    all_exprds_fcc
    
  }

  
  outFile <- file.path(outFolder, "all_permut_fcc.Rdata")
  save(all_permut_fcc, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_permut_fcc.Rdata")
  all_permut_fcc <- get(load(outFile))
}


outFile <- file.path(outFolder, paste0("allDS_cmp_ratioDown_meanCorrPermut.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot_multiDens(
 list( ratioDown_meanRL_v2 = unlist(lapply(all_permut_fcc, function(sub) lapply(sub, function(x) x[["ratioDown_meanRL"]]))),
  ratioDown_RorL_v2_either = unlist(lapply(all_permut_fcc, function(sub) lapply(sub, function(x) x[["ratioDown_RorL"]]))),
#  ratioDown_meanRL_withT_v3 = unlist(lapply(all_permut_fcc, function(sub) lapply(sub, function(x) x[["ratioDown_meanRL_withTAD"]]))),
#  ratioDown_RorL_withT_v3_either = unlist(lapply(all_permut_fcc, function(sub) lapply(sub, function(x) x[["ratioDown_RorL_withTAD"]]))),
  ratioDown_RorL_replace_v4_either = unlist(lapply(all_permut_fcc, function(sub) lapply(sub, function(x) x[["ratioDown_RorL_replace"]])))
 ), legPos="topleft"
)
cat(paste0("... written: ", outFile, "\n"))
