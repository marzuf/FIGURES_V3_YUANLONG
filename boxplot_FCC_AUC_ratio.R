
# Rscript boxplot_FCC_AUC_ratio.R


plotType <-  "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

outFolder <- file.path("BOXPLOT_FCC_AUC_RATIO")
dir.create(outFolder, recursive = TRUE)

obs_dt <- get(load("BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"))
rd_midpos_dt <- get(load("BARPLOT_FCC_AUC_RATIO_RANDOMMIDPOS//all_dt.Rdata"))
rd_midposdisc_dt <- get(load("BARPLOT_FCC_AUC_RATIO_RANDOMMIDPOSDISC///all_dt.Rdata"))
rd_midposstrict_dt <- get(load("BARPLOT_FCC_AUC_RATIO_RANDOMMIDPOSSTRICT///all_dt.Rdata"))

all_dt <- do.call(rbind, list(obs_dt, rd_midposdisc_dt, rd_midpos_dt, rd_midposstrict_dt))

all_dt$hicds_lab <- gsub(".+_(.+?)_40kb","\\1",  all_dt$hicds)
all_dt$hicds_lab <- ifelse(grepl("RANDOM",all_dt$hicds_lab) | grepl("PERMUT", all_dt$hicds_lab),all_dt$hicds_lab,  "OBSERVED" )

stopifnot(diff(table(all_dt$hicds_lab)) == 0)

nDS <- as.numeric(unique(table(all_dt$hicds_lab)))

outFile <- file.path(outFolder, paste0("fcc_auc_ratio_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
par(mar = par()$mar + c(9,0,0,0))
boxplot(fcc_auc~hicds_lab, data = all_dt, main = "FCC AUC ratio", xlab="", ylab="FCC AUC ratio", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
mtext(side=3, text = paste0("all datasets - n=", nDS))
foo <- dev.off()




outFile <- file.path(outFolder, paste0("coexpr_auc_ratio_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
par(mar = par()$mar + c(9,0,0,0))
boxplot(coexpr_auc~hicds_lab, data = all_dt, main = "Coexpr. AUC ratio", xlab="", ylab="Coexpr. AUC ratio", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
mtext(side=3, text = paste0("all datasets - n=", nDS))
foo <- dev.off()


#load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/AUC_COEXPRDIST_WITHFAM_SORTNODUP/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas_hgnc/hgnc_family_short/auc_values.Rdata")
#obs_dt$coexpr_auc[1] == auc_values[["auc_ratio_same_over_diff_distVect"]]
##[1] TRUE

