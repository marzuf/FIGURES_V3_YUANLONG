
# Rscript scatterplot_fcc_fract_withRandom_either.R

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)

require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

myWidth <- myWidth * 1.2

outFolder <- "SCATTERPLOT_FCC_FRACT_WITHRANDOM_EITHER"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

fcc_fract <- seq(from=-1, to=1, by=0.25)
# fcc_fract_names <- paste0("FCC > ", fcc_fract[1:(length(fcc_fract)-1)], " and FCC <= ",fcc_fract[2:length(fcc_fract)])
fcc_fract_names <- paste0("FCC \u2208 ]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
fcc_fract_names <- paste0("]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")

fcc_fract_names[fcc_fract_names == "]-1, -0.75]"] <- "[-1, -0.75]"


# fract_sort <- "FCC > 0.75 and FCC <= 1"
fract_sort <- fcc_fract_names[length(fcc_fract_names)]

ggsci_pal <- "lancet"
ggsci_subpal <- ""

legTitle <- "FCC ranges:"
fractBarSubTitle <- "AUC ratios:\n"
fractBarTitle <- "Fold-change concordance scores"

auc_ratio_file <- file.path("BARPLOT_FCC_AUC_RATIO", "all_dt.Rdata")
stopifnot(file.exists(auc_ratio_file))

rd_fcc_folder <- file.path("RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V2_EITHER")

# SETTINGS FOR THE PLOTS
dotPch <- 19
dotCex <- 0.7
mainCex <- 1.2
rd_opc <- 0.3


if(buildData){
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    cat(paste0("... start: ", hicds, "\n"))
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
      fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      all_fcc <- get(load(fcc_file))

      fcc_hist <- hist(all_fcc, breaks=fcc_fract)$counts
      names(fcc_hist) <- fcc_fract_names
      
      fcc_hist <- fcc_hist/length(all_fcc)
      stopifnot(sum(fcc_hist) == 1)
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        intervalFCC = names(fcc_hist),
        countFCC = as.numeric(fcc_hist),
        stringsAsFactors = FALSE
      )
    }
    exprds_dt
  }
  outFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  auc_fract_file <- outFile
  save(all_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  inFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  auc_fract_file <- inFile
  all_dt <- get(load(inFile))
  # load("SCATTERPLOT_FCC_FRACT_WITHRANDOM/all_dt.Rdata")
}

######################################################################################
# FRACT and FCC AUC - discrete x-axis
######################################################################################


# auc_fract_dt <- get(load("SCATTERPLOT_FCC_FRACT_WITHRANDOM/all_dt.Rdata"))
# auc_ratio_dt <- get(load("BARPLOT_FCC_AUC_RATIO//all_dt.Rdata"))

auc_fract_dt <- get(load(auc_fract_file))
auc_ratio_dt <- get(load(auc_ratio_file))

auc_ratio_dt$dataset <- paste0(auc_ratio_dt$hicds, "\n", auc_ratio_dt$exprds)
auc_ratio_dt <- auc_ratio_dt[order(auc_ratio_dt$fcc_auc, decreasing = TRUE),]
auc_ds_order <- auc_ratio_dt$dataset
auc_ds_cols <- all_cols[all_cmps[paste0(auc_ratio_dt$exprds)]]

auc_fract_dt$dataset <- paste0(auc_fract_dt$hicds, "\n", auc_fract_dt$exprds)

auc_fract_ratio_dt <- merge(auc_fract_dt, auc_ratio_dt, by=c("hicds", "exprds", "dataset"), all.x=TRUE, all.y=TRUE)
stopifnot(!is.na(auc_fract_ratio_dt))
auc_fract_ratio_dt$dataset <- factor(auc_fract_ratio_dt$dataset, levels=auc_ds_order)
stopifnot(!is.na(auc_fract_ratio_dt$dataset))

plot_dt <- auc_fract_ratio_dt
plot_dt$rank <- as.numeric(plot_dt$dataset)

plot_dt$intervalFCC <- factor(plot_dt$intervalFCC, levels = fcc_fract_names)
stopifnot(!is.na(plot_dt$intervalFCC))

maxRank <- max(plot_dt$rank) 
nDS <- length(unique(auc_ds_order))

myPals <-  eval(parse(text=paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(unique(fcc_fract_names)))
myPals <- rev(myPals)

my_xlab <- "Datasets (ranked by decreasing FCC AUC ratio)"

my_ylab <- "Fract. of TADs"


#all_rd_types <- c("", "_meanRL")
all_rd_types <- c("_RorL")
rd_type=all_rd_types[1]

for(rd_type in all_rd_types) {
  
  if(rd_type == "") {
    plotTit <- paste0("all datasets (n=", nDS, ") (RandL)")  
  } else {
    plotTit <- paste0("all datasets (n=", nDS, ") (", gsub("_", "", rd_type),")")  
  }
  
  
  
  rdFile <- file.path(rd_fcc_folder, paste0("rd", rd_type, "_lm_intervalFCC.Rdata"))
  # load("RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT/rd_lm_intervalFCC.Rdata")
  rd_lm_intervalFCC <- get(load(rdFile))
  
  rd_y_max <- max(unlist(lapply(rd_lm_intervalFCC, function(x) {
    x$coefficients[1] + 1:maxRank * x$coefficients[2]
  })))
  rd_y_min <- min(unlist(lapply(rd_lm_intervalFCC, function(x) {
    x$coefficients[1] + 1:maxRank * x$coefficients[2]
  })))
  y_max <- max(c(plot_dt$countFCC, rd_y_max))
  y_min <- min(c(plot_dt$countFCC, rd_y_min))

  
  ###################### PLOT WITH THE YLABS 
  
  outFile <- file.path(outFolder, paste0("rd", rd_type, "_fcc_fract_withRandom_withLabs_scatterplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  dev.control(displaylist="enable")
  par(bty="L", family = fontFamily)
  
  plot(NULL,
       xlim=c(1, maxRank),
       ylim = c(y_min, y_max),
       xlab = paste0(my_xlab),
       ylab = paste0(my_ylab),
       axes = FALSE,
       main = plotTit,
       cex.main = mainCex,
       cex.lab = labCex,
       cex.axis = axisCex
       )
  box(bty="L")
  
  axis(2, at = seq(0, y_max + 0.05, by=0.05) , lwd=0, lwd.ticks=1)
  axis(1, at = seq(1, maxRank, by=1) , lwd=0, labels=F,lwd.ticks=1, tck=-0.01)
  
  # axis(1, at = 1:maxRank,
  #     labels = rep("", maxRank)
  #      )
  # axis(1, 
  #      labels = FALSE,
  #      lwd.ticks = -1
  # )
  
  
  i_fract = 8
  for(i_fract in 1:length(fcc_fract_names)) {
    
    dotCol <- myPals[i_fract]
    sub_plot_dt <- plot_dt[as.numeric(plot_dt$intervalFCC) == i_fract,]
    
    points(x = sub_plot_dt$rank, y = sub_plot_dt$countFCC,
           pch = dotPch,
           cex = dotCex,
           col = dotCol)
    
    abline(
      lm(sub_plot_dt$countFCC ~ sub_plot_dt$rank), 
      lty = 2,
      col = dotCol
    )
    
    # load and plot the random model
    rd_lm <- rd_lm_intervalFCC[[paste0(fcc_fract_names[i_fract])]]
    
    # Y=255 - P*(255-X) where X is a RGB number, P=opacity (0...1), Y=new RGB number 
    rd_col <- (255 - rd_opc*(255-c(col2rgb(dotCol))))/255
    
    # curve(rd_lm$coefficients[1] * x + rd_lm$coefficients[2], add = TRUE, col = dotCol, lty=1)
    curve(rd_lm$coefficients[1] + x*rd_lm$coefficients[2], add = TRUE, col = rgb(matrix(rd_col, ncol=3)), lty=1, lwd=5)

    
  } # end-for iterating FCC fract

  # blackcol_leg <- (255 - 0.7*(255-c(col2rgb("black"))))/255
  blackcol_leg <- "black"
  # legend(
  #   # "topleft",
  #   "top",
  #   # inset = c(-0.05, -0.05),
  #   inset = c(-0.0, -0.06),
  #   xpd = TRUE,
  #   # horiz = TRUE,
  #   pch=dotPch,
  #   ncol = length(fcc_fract_names)/2,
  #   legend=c(rev(levels(plot_dt$intervalFCC))),
  #   lty=1,
  #   lwd = 2,
  #   col = rev(myPals),
  #   bty="n"
  # )
  legend(
    # "topleft",
    "top",
    # inset = c(-0.05, -0.05),
    inset = c(-0.0, -0.08),
    xpd = TRUE,
    # horiz = TRUE,
    pch=c(rep(dotPch, length(fcc_fract_names)), -1, -1),
    ncol = (length(fcc_fract_names)/2+1),
    legend=c(rev(levels(plot_dt$intervalFCC)), "fit permut.", "fit obs."),
    lty=c(rep(1, length(fcc_fract_names)), 1, 2),
    lwd = c(rep(2, length(fcc_fract_names)), 5, 1),
    col = c(rev(myPals), blackcol_leg, blackcol_leg),
    bty="n"
  )
  scatP <- recordPlot()

  # mtext(1, text = auc_ds_order, col = auc_ds_cols, at=1:maxRank, las=2, cex=0.4, line=0.8)
  mtext(1, text = auc_ds_order, col = auc_ds_cols, at=1:maxRank, las=2, cex=0.4, adj=0, line=0.8)
  
  
  invisible(dev.off())
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0("rd", rd_type, "_fcc_fract_withRandom_withSymb_scatterplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  replayPlot(scatP) 
  mtext(1, text = rep(labsymbol, maxRank), col = auc_ds_cols, at=1:maxRank, las=2, cex=0.7, line=0.8)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
} # end-for iterating rd_type



