
options(scipen=100)

# Rscript random_FCC_AUC_ratio_onePerm.R
script_name <- "random_FCC_AUC_ratio_onePerm.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildData <- FALSE

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

outFolder <- file.path("RANDOM_FCC_AUC_RATIO_ONEPERM")
dir.create(outFolder, recursive = TRUE)

runFolder <- file.path("..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA")

permutFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")

plotMargin <- c(1,2,1,1)
legTitle <- "FCC ranges:"
fractBarSubTitle <- "AUC ratios:\n"
fractBarTitle_main <- "Fold-change concordance scores - PERMUT DATA"

ggsci_pal <- "lancet"
ggsci_subpal <- ""


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "svg"
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")


script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script8_name <- "8cOnlyFCC_runAllDown"
permut_script <- "5sameNbr_runPermutationsCorr"

fcc_fract <- seq(from=-1, to=1, by=0.25)
# fcc_fract_names <- paste0("FCC > ", fcc_fract[1:(length(fcc_fract)-1)], " and FCC <= ",fcc_fract[2:length(fcc_fract)])
fcc_fract_names <- paste0("FCC \u2208 ]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
fcc_fract_names <- paste0("]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
fcc_fract_names[fcc_fract_names == "]-1, -0.75]"] <- "[-1, -0.75]"
# fract_sort <- "FCC > 0.75 and FCC <= 1"
fract_sort <- fcc_fract_names[length(fcc_fract_names)]


all_hicds <- list.files(file.path(permutFolder))
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(permutFolder, x)))


get_fcc <- function(fc_vect) {
  # (2* sum(all_tad_fc < 0)/length(all_tad_fc) -1) *  (2* sum(abs(all_tad_fc[all_tad_fc<0]))/sum(abs(all_tad_fc)) -1)
  # diffRatioFC <- (sum(abs(all_tad_fc[all_tad_fc < 0])) - sum(abs(all_tad_fc[all_tad_fc > 0])) ) /sum(abs(all_tad_fc))
  # ratioDown <- sum(all_tad_fc < 0)/length(all_tad_fc)
  # stopifnot(ratioDown >=0 & ratioDown <=1)
  # stopifnot(diffRatioFC >=-1 & diffRatioFC <=1)
  # prodRatio <- 2*(ratioDown - 0.5) * diffRatioFC
  # stopifnot(all(prodRatio >= -1 & prodRatio <= 1))
  (2* sum(fc_vect < 0)/length(fc_vect) -1) *  (2* sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect)) -1)
}

get_ratioDown <- function(fc_vect) {
   sum(fc_vect < 0)/length(fc_vect) 
}

get_ratioFC <- function(fc_vect) {
  # (2* sum(all_tad_fc < 0)/length(all_tad_fc) -1) *  (2* sum(abs(all_tad_fc[all_tad_fc<0]))/sum(abs(all_tad_fc)) -1)
  # diffRatioFC <- (sum(abs(all_tad_fc[all_tad_fc < 0])) - sum(abs(all_tad_fc[all_tad_fc > 0])) ) /sum(abs(all_tad_fc))
  # ratioDown <- sum(all_tad_fc < 0)/length(all_tad_fc)
  # stopifnot(ratioDown >=0 & ratioDown <=1)
  # stopifnot(diffRatioFC >=-1 & diffRatioFC <=1)
  # prodRatio <- 2*(ratioDown - 0.5) * diffRatioFC
  # stopifnot(all(prodRatio >= -1 & prodRatio <= 1))
  sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect))
}


set.seed(742020)

if(buildData) {
  all_permut_fcc <- foreach(hicds = all_hicds) %do%{
    all_exprds_fcc <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      permut_data <- get(load(file.path(permutFolder, hicds, exprds, script8_name, "prodSignedRatio_permDT.Rdata")))
      # keep only 1 permut
      i_keep <- sample(1:ncol(permut_data), size=1)
      fcc_onePermut <- permut_data[,i_keep]
      
      
      names(fcc_onePermut) <- rownames(permut_data)

            
      
      # compute the FCC AUC for the obs. data

      obs_fcc <- get(load(file.path(permutFolder, hicds, exprds, script8_name, "all_obs_prodSignedRatio.Rdata")))
      obs_fcc_sorted <- sort(obs_fcc, decreasing = TRUE)
      
      
      rd_fcc_onePermut <- fcc_onePermut
      stopifnot(length(rd_fcc_onePermut) == length(obs_fcc_sorted))
      rd_fcc_onePermut <- na.omit(rd_fcc_onePermut)
      
      rd_fcc_onePermut_hist <- hist(rd_fcc_onePermut, breaks=fcc_fract)$counts
      names(rd_fcc_onePermut_hist) <- fcc_fract_names
      rd_fcc_onePermut_ratio_hist <- rd_fcc_onePermut_hist/length(rd_fcc_onePermut)
      stopifnot(abs(sum(rd_fcc_onePermut_ratio_hist) - 1) < 10^-4)
      
      rd_fcc_onePermut_sorted <- sort(rd_fcc_onePermut, decreasing = TRUE)
      
      nTot_onePermut <- min(c(length(rd_fcc_onePermut_sorted), length(obs_fcc_sorted)))
      x_val_onePermut <- 1:nTot_onePermut
      
      
      cumsum_obs <- cumsum(abs(obs_fcc_sorted[1:nTot_onePermut])) # updated here 20.03.2020 ! should be abs !
      cumsum_rd <- cumsum(abs(rd_fcc_onePermut_sorted[1:nTot_onePermut]))
      
      auc_obs_onePermut <- auc(x = x_val_onePermut, y = cumsum_obs)
      auc_rd_onePermut <- auc(x = x_val_onePermut, y = cumsum_rd)
      
      auc_ratio_rd_onePermut <- auc_obs_onePermut/auc_rd_onePermut
      
      

      
      list(
        auc_ratio_rd_onePermut = auc_ratio_rd_onePermut,
        rd_fcc_onePermut_ratio_hist = rd_fcc_onePermut_ratio_hist
      )
    } # end-foreach exprds
    names(all_exprds_fcc) <- all_exprds[[paste0(hicds)]]
    all_exprds_fcc
  } # end-foreach hicds
  names(all_permut_fcc) <- all_hicds 
  
  outFile <- file.path(outFolder, "all_permut_fcc.Rdata")
  save(all_permut_fcc, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  # load("RANDOM_FCC_AUC_RATIO_ONEPERM/all_permut_fcc.Rdata")
  outFile <- file.path(outFolder, "all_permut_fcc.Rdata")
  all_permut_fcc <- get(load(outFile))
}

obs_auc_dt <- get(load("BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"))


all_permut_fcc_ul <- unlist(all_permut_fcc, recursive = FALSE)

# for each dataset -> 
all_rd_types <- c("_onePermut")
rd_type=all_rd_types[1]

for(rd_type in all_rd_types) {
  
  if(rd_type=="") {
    fractBarTitle <- paste0(fractBarTitle_main, " (RandL)")
    scatterTit <- paste0("PERMUT DATA (RandL)")
    cmpTit <- "(RandL)"
  }else {
    fractBarTitle <- paste0(fractBarTitle_main, " (", gsub("_", "", rd_type), ")")  
    scatterTit <- paste0("PERMG2T (1 permut)")
    cmpTit <- paste0("(", gsub("_", "", rd_type), ")")
  }
  
  
  rd_auc_fcc_ratio_dt <- data.frame(do.call(rbind, lapply(all_permut_fcc_ul, function(x) x[[paste0("auc_ratio_rd", rd_type)]])))
  colnames(rd_auc_fcc_ratio_dt)[1] <- "rd_fcc_auc"
  rd_auc_fcc_ratio_dt$hicds <- gsub("^(.+)\\.TCGA.+$", "\\1", rownames(rd_auc_fcc_ratio_dt))
  rd_auc_fcc_ratio_dt$exprds <- gsub("^.+\\.(TCGA.+)$", "\\1", rownames(rd_auc_fcc_ratio_dt))
  rownames(rd_auc_fcc_ratio_dt) <- NULL
  rd_auc_fcc_ratio_dt <- rd_auc_fcc_ratio_dt[order(rd_auc_fcc_ratio_dt$rd_fcc_auc, decreasing = TRUE),]
  rd_auc_fcc_ratio_dt$dataset <- paste0(rd_auc_fcc_ratio_dt$hicds, "\n", rd_auc_fcc_ratio_dt$exprds)
  
  mycols <- all_cols[all_cmps[paste0(rd_auc_fcc_ratio_dt$exprds)]]
  
  nDS <- length(unique(rd_auc_fcc_ratio_dt$dataset))
  myds <- as.character(unique(file.path(rd_auc_fcc_ratio_dt$hicds, rd_auc_fcc_ratio_dt$exprds)))
  countCmp <- setNames(as.numeric(table(all_cmps[basename(myds) ])), names(table(all_cmps[basename(myds) ])))
  
  legDT <- data.frame(cmpType = names(all_cols), cmpColor = as.character(all_cols))
  legDT <- legDT[rev(order(legDT$cmpType)),]
  legDT$count <- countCmp[legDT$cmpType]
  legDT$legLabel <- paste0(legDT$cmpType, " (", legDT$count, ")")
  
  
  outFile <- file.path(outFolder, paste0("rd", rd_type, "_auc_fcc_ratio_dt.Rdata"))
  save(rd_auc_fcc_ratio_dt, file = outFile)
  
  rd_auc_fcc_ratio_dt$dataset <- paste0("rd_", rd_auc_fcc_ratio_dt$dataset)
  auc_ds_order <- as.character(rd_auc_fcc_ratio_dt$dataset)
  
  
  rd_auc_fcc_ratio_fract_dt <- data.frame(do.call(rbind, lapply(all_permut_fcc_ul, function(x) x[[paste0("rd_fcc", rd_type, "_ratio_hist")]])), check.names = FALSE)
  rd_auc_fcc_ratio_fract_dt$hicds <- gsub("^(.+)\\.TCGA.+$", "\\1", rownames(rd_auc_fcc_ratio_fract_dt))
  rd_auc_fcc_ratio_fract_dt$exprds <- gsub("^.+\\.(TCGA.+)$", "\\1", rownames(rd_auc_fcc_ratio_fract_dt))
  rownames(rd_auc_fcc_ratio_fract_dt) <- NULL
  rd_auc_fcc_ratio_fract_dt$dataset <- paste0(rd_auc_fcc_ratio_fract_dt$hicds, "\n", rd_auc_fcc_ratio_fract_dt$exprds)
  rd_auc_fcc_ratio_fract_dt <- melt(rd_auc_fcc_ratio_fract_dt, id =c("hicds", "exprds", "dataset"))
  colnames(rd_auc_fcc_ratio_fract_dt)[colnames(rd_auc_fcc_ratio_fract_dt) == "variable"] <- "rd_intervalFCC"  
  colnames(rd_auc_fcc_ratio_fract_dt)[colnames(rd_auc_fcc_ratio_fract_dt) == "value"] <- "rd_countFCC"  
  
  rd_auc_fcc_ratio_fract_dt$dataset <- paste0("rd_", rd_auc_fcc_ratio_fract_dt$dataset)
  rd_auc_fcc_ratio_fract_dt$dataset <- factor(rd_auc_fcc_ratio_fract_dt$dataset, levels=auc_ds_order)
  stopifnot(!is.na(rd_auc_fcc_ratio_fract_dt$dataset))
  
  
  rd_auc_fcc_ratio_fract_dt$rd_intervalFCC <- factor(rd_auc_fcc_ratio_fract_dt$rd_intervalFCC, levels = rev(fcc_fract_names))
  stopifnot(!is.na(rd_auc_fcc_ratio_fract_dt$rd_intervalFCC))
  
  ############################################################### 
  ############################################################### cmp FCC AUC ranking
  ############################################################### 
  
  cmp_dt <- merge(rd_auc_fcc_ratio_dt, obs_auc_dt, by=c("hicds", "exprds"))
  
  cmp_dt$dotCols <- all_cols[all_cmps[paste0(cmp_dt$exprds)]]
  
  outFile <- file.path(outFolder, paste0("rd", rd_type, "_rd_vs_obs_scatterplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="l")
  plot(
    x = cmp_dt$fcc_auc,
    y = cmp_dt$rd_fcc_auc,
    pch = 16,
    cex = 1,
    main = "Permut. vs. obs. FCC AUC ratios",
    xlab = "Obs. FCC AUC ratio",
    ylab = "Permut FCC AUC ratio",
    col = cmp_dt$dotCols
  )
  mtext(side = 3, text = cmpTit)
  legend(
    "topright",
    legend = legDT$legLabel,
    col = as.character(legDT$cmpColor),
    pch=16,
    bty="n"
  )
  addCorr(legPos = "bottomleft", x= cmp_dt$fcc_auc, y = cmp_dt$rd_fcc_auc, bty="n")
  abline(lm(cmp_dt$rd_fcc_auc~cmp_dt$fcc_auc), lty=2, col="grey")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  ############################################################### 
  ############################################################### BARPLOT
  ############################################################### 
  
  
  rd_fract_plot_with_lab <- ggplot(rd_auc_fcc_ratio_fract_dt, aes(x=dataset, y=rd_countFCC, fill=rd_intervalFCC, color=rd_intervalFCC)) + 
    geom_bar(position="stack", stat="identity") +
    ggtitle(paste0(fractBarTitle), 
            subtitle = paste0(fractBarSubTitle) )+
    # subtitle = "(all datasets)")+
    scale_x_discrete(name=paste0("(all datasets - n=", nDS, ")"))+
    labs(fill=paste0(legTitle)) + 
    guides(color=FALSE)+
    # coord_cartesian(expand=FALSE)+
    coord_cartesian(clip = 'off', expand=FALSE) +
    scale_y_continuous(name=paste0("Fraction of TADs"),
                       limits = c(0,1), 
                       breaks = seq(from=0, to=1, by=0.1),
                       labels = seq(from=0, to=1, by=0.1))+
    eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
    eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
    theme( # Increase size of axis lines
      plot.margin = unit(plotMargin, "lines"), # top, right, bottom, and left 
      plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
      plot.subtitle = element_text(hjust = 0, vjust=2,face = "italic", size =8, family=fontFamily, lineheight = 1.75),
      panel.grid = element_blank(),
      # panel.grid.major.y = element_line(colour = "grey"),
      # panel.grid.minor.y = element_line(colour = "grey"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .2, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
      axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90, family=fontFamily),
      axis.title.y = element_text(color="black", size=14, family=fontFamily),
      axis.title.x = element_text(color="black", size=14, family=fontFamily),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank(),
      legend.title = element_text(face="bold")
    ) +
    geom_text(data=rd_auc_fcc_ratio_dt, aes(x = rd_auc_fcc_ratio_dt$dataset, y=1, 
                          label=sprintf("%.2f", rd_auc_fcc_ratio_dt$rd_fcc_auc)),
              inherit.aes=FALSE, angle=90, size=3, 
              vjust=0.5, hjust=0)+
    theme(
      legend.justification = c("right", "top")
    )+ 
    geom_text(data=legDT, aes(label = legDT$legLabel, x = 59, y =c(0, 0.05, 0.1)),
              vjust = 0, hjust=0,
              inherit.aes = FALSE, color = legDT$cmpColor)
  
  
  
  outFile <- file.path(outFolder, paste0("rd", rd_type, "_all_ds_fcc_fract_scores_withLabs_barplot.", plotType))
  ggsave(plot = rd_fract_plot_with_lab, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*2)
  cat(paste0("... written: ", outFile, "\n"))
  
  ############################################################### 
  ############################################################### SCATTERPLOT
  ############################################################### 
  
  rd_auc_fcc_ratio_dt <- data.frame(do.call(rbind, lapply(all_permut_fcc_ul, function(x) x[[paste0("auc_ratio_rd", rd_type)]])))
  colnames(rd_auc_fcc_ratio_dt)[1] <- "rd_fcc_auc"
  rd_auc_fcc_ratio_dt$hicds <- gsub("^(.+)\\.TCGA.+$", "\\1", rownames(rd_auc_fcc_ratio_dt))
  rd_auc_fcc_ratio_dt$exprds <- gsub("^.+\\.(TCGA.+)$", "\\1", rownames(rd_auc_fcc_ratio_dt))
  rownames(rd_auc_fcc_ratio_dt) <- NULL
  rd_auc_fcc_ratio_dt <- rd_auc_fcc_ratio_dt[order(rd_auc_fcc_ratio_dt$rd_fcc_auc, decreasing = TRUE),]
  rd_auc_fcc_ratio_dt$dataset <- paste0(rd_auc_fcc_ratio_dt$hicds, "\n", rd_auc_fcc_ratio_dt$exprds)
  
  
  
  rd_auc_fcc_ratio_fract_dt <- data.frame(do.call(rbind, lapply(all_permut_fcc_ul, function(x) x[[paste0("rd_fcc", rd_type, "_ratio_hist")]])), check.names = FALSE)
  rd_auc_fcc_ratio_fract_dt$hicds <- gsub("^(.+)\\.TCGA.+$", "\\1", rownames(rd_auc_fcc_ratio_fract_dt))
  rd_auc_fcc_ratio_fract_dt$exprds <- gsub("^.+\\.(TCGA.+)$", "\\1", rownames(rd_auc_fcc_ratio_fract_dt))
  rownames(rd_auc_fcc_ratio_fract_dt) <- NULL
  rd_auc_fcc_ratio_fract_dt$dataset <- paste0(rd_auc_fcc_ratio_fract_dt$hicds, "\n", rd_auc_fcc_ratio_fract_dt$exprds)
  rd_auc_fcc_ratio_fract_dt <- melt(rd_auc_fcc_ratio_fract_dt, id =c("hicds", "exprds", "dataset"))
  colnames(rd_auc_fcc_ratio_fract_dt)[colnames(rd_auc_fcc_ratio_fract_dt) == "variable"] <- "rd_intervalFCC"  
  colnames(rd_auc_fcc_ratio_fract_dt)[colnames(rd_auc_fcc_ratio_fract_dt) == "value"] <- "rd_countFCC"  
  
  rd_auc_fcc_ratio_fract_dt$dataset <- paste0("rd_", rd_auc_fcc_ratio_fract_dt$dataset)
  rd_auc_fcc_ratio_fract_dt$dataset <- factor(rd_auc_fcc_ratio_fract_dt$dataset, levels=auc_ds_order)
  stopifnot(!is.na(rd_auc_fcc_ratio_fract_dt$dataset))
  
  rd_auc_fcc_ratio_fract_dt$ds_rank <- as.numeric(rd_auc_fcc_ratio_fract_dt$dataset)
  
    
  mycols_scat <- all_cols[all_cmps[paste0(rd_auc_fcc_ratio_dt$exprds)]]
  xax_exp <- 0.01
  myPals <-  eval(parse(text=paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(unique(rd_auc_fcc_ratio_fract_dt$rd_intervalFCC)))
  xbreaks <- sort(unique(rd_auc_fcc_ratio_fract_dt$ds_rank))
  
  rd_auc_fcc_ratio_fract_dt$rd_intervalFCC <- factor(rd_auc_fcc_ratio_fract_dt$rd_intervalFCC, levels = rev(fcc_fract_names))
  stopifnot(!is.na(rd_auc_fcc_ratio_fract_dt$rd_intervalFCC))
  
  
  rd_scatPlot <- ggscatter(rd_auc_fcc_ratio_fract_dt, 
                        title = paste0(scatterTit, " (n=", length(unique(file.path(rd_auc_fcc_ratio_fract_dt$hicds, rd_auc_fcc_ratio_fract_dt$exprds))), ")"),
                        x = "ds_rank", 
                        y = "rd_countFCC",
                        color = "rd_intervalFCC",
                        xlab = "Datasets (ranked by FCC AUC ratio)",
                        ylab = "Ratio of TADs",
                        palette = myPals) +
    labs(color=legTitle)+
    geom_smooth(aes(color = rd_intervalFCC),method = "lm", linetype=2, se=F)+
    
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_x_continuous(labels = auc_ds_order, breaks=xbreaks, expand=c(xax_exp,xax_exp))+
    
    theme( # Increase size of axis lines
      plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14, family=fontFamily),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
      axis.text.x = element_text(color=mycols_scat, hjust=1,vjust = 0.5, size=7, family=fontFamily, angle=90),
      axis.title.y = element_text(color="black", size=14, family=fontFamily),
      axis.title.x = element_text(color="black", size=14, family=fontFamily),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.title = element_text(face="bold", family=fontFamily)
    )
  
  outFile <- file.path(outFolder, paste0("rd", rd_type, "_all_ds_fcc_fract_scores_nbrSignifs_auc_scatterplot_xcont_withLabs.", plotType))
  ggsave(plot = rd_scatPlot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  # save the lm curves to be able to re-draw the curves
  # -> for each of the ranges, fit lm
  
  
  rd_fccRange = fcc_fract_names[1]
  rd_lm_intervalFCC <- foreach(rd_fccRange = fcc_fract_names) %dopar% {
    sub_dt <- rd_auc_fcc_ratio_fract_dt[as.character(rd_auc_fcc_ratio_fract_dt$rd_intervalFCC) == rd_fccRange,]
    lm(sub_dt$rd_countFCC ~ sub_dt$ds_rank)
  }
  names(rd_lm_intervalFCC) <- fcc_fract_names
  outFile <- file.path(outFolder, paste0("rd", rd_type, "_lm_intervalFCC.Rdata"))
  save(rd_lm_intervalFCC, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
}


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

