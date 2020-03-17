
# Rscript barplot_with_FCC_fract_randommidposstrict.R

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)

require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")


outFolder <- "BARPLOT_WITH_FCC_FRACT_RANDOMMIDPOSSTRICT"
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

plotMargin <- c(1,2,1,1)

obs_ratio_file <- "BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"
stopifnot(file.exists(obs_ratio_file))
rd_ratio_file <- "BARPLOT_FCC_AUC_RATIO_RANDOMMIDPOSSTRICT/all_dt.Rdata"
stopifnot(file.exists(rd_ratio_file))

auc_ratio_file <- file.path("BARPLOT_FCC_AUC_RATIO_RANDOMMIDPOSSTRICT", "all_dt.Rdata")
stopifnot(file.exists(auc_ratio_file))
x <- get(load(auc_ratio_file))
x$dataset <- paste0(x$hicds, "\n", x$exprds)
x <- x[order(x$fcc_auc, decreasing=TRUE),]
fcc_ds_order <- x$dataset
fcc_auc_sort <- TRUE

signif_file <- file.path(runFolder, "CREATE_FINAL_TABLE_RANDOM/all_result_dt.Rdata")
stopifnot(file.exists(signif_file))

all_hicds <- all_hicds[grepl("_RANDOMMIDPOSSTRICT_", all_hicds)]

if(buildData){
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    cat(paste0("... start: ", hicds, "\n"))
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
      fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      all_fcc <- get(load(fcc_file))
      
      # [1] -1.00 -0.75   -0.50       -0.25   0.00    0.25  0.50  0.75  1.00
      # [1]   0   0        10        411      571      355  97 350
      # > sum(all_fcc > 0.75 & all_fcc <= 1)
      # [1] 350
      # > sum(all_fcc > 0.5 & all_fcc <= 0.75)
      # [1] 97
      # sum(all_fcc > -0.5 & all_fcc <= -0.25)
      # [1] 10
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
  # load("BARPLOT_WITH_FCC_FRACT/all_dt.Rdata")
}

all_dt$intervalFCC <- factor(all_dt$intervalFCC, levels = rev(fcc_fract_names))
stopifnot(!is.na(all_dt$intervalFCC))

all_dt$dataset <- paste0(all_dt$hicds, "\n", all_dt$exprds)
all_dt <- all_dt[order(all_dt$countFCC, decreasing = TRUE),]
ds_order <- all_dt$dataset[all_dt$intervalFCC == fract_sort]
exprds_order <- all_dt$exprds[all_dt$intervalFCC == fract_sort]
exprds_order_fract <- exprds_order

all_dt$labSymb <- labsymbol

if(fcc_auc_sort) {
  all_dt$dataset <- factor(all_dt$dataset, levels=fcc_ds_order)
  stopifnot(!is.na(all_dt$dataset))
  
} else {
  all_dt$dataset <- factor(all_dt$dataset, levels=ds_order)
  stopifnot(!is.na(all_dt$dataset))
}
tmp1 <- all_dt
tmp1 <- tmp1[order(as.numeric(tmp1$dataset)),]
exprds_order <- as.character(tmp1$exprds[tmp1$intervalFCC == fract_sort])

nDS <- length(unique(all_dt$dataset))

myds <- as.character(unique(file.path(all_dt$hicds, all_dt$exprds)))
countCmp <- setNames(as.numeric(table(all_cmps[basename(myds) ])), names(table(all_cmps[basename(myds) ])))

legDT <- data.frame(cmpType = names(all_cols), cmpColor = as.character(all_cols))
legDT <- legDT[rev(order(legDT$cmpType)),]
legDT$count <- countCmp[legDT$cmpType]
legDT$legLabel <- paste0(legDT$cmpType, " (", legDT$count, ")")

mycols <- all_cols[all_cmps[exprds_order]]
# tmp_dt <- aggregate(countFCC~ intervalFCC, data=all_dt, FUN=sum) # none has zero

mycols_fract <- all_cols[all_cmps[exprds_order_fract]]

fract_plot_tmp <- ggplot(all_dt, aes(x=dataset, y=countFCC, fill=intervalFCC, color=intervalFCC)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle(paste0(fractBarTitle), 
          subtitle = paste0(fractBarSubTitle) )+
  labs(fill=paste0(legTitle)) + 
  guides(color=FALSE)+
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
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  ) +
  geom_text(data=x, aes(x = x$dataset, y=1, 
                        label=sprintf("%.2f", x$fcc_auc)),
            inherit.aes=FALSE, angle=90, size=3, 
            vjust=0.5, hjust=0)+
  theme(
    legend.justification = c("right", "top")
  )+ 
  geom_text(data=legDT, aes(label = legDT$legLabel, x = 59, y =c(0, 0.05, 0.1)),
            vjust = 0, hjust=0,
            inherit.aes = FALSE, color = legDT$cmpColor)

fract_plot_with_lab <- fract_plot_tmp +
  scale_x_discrete(name=paste0("(all datasets - n=", nDS, ")"))+
  theme(axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90, family=fontFamily))

outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withLabs_barplot.", plotType))
ggsave(plot = fract_plot_with_lab, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))


fract_plot_with_symb <- fract_plot_tmp + 
scale_x_discrete(labels=all_dt$labSymb, name=paste0("(all datasets - n=", nDS, ")"))+
  theme(axis.text.x = element_text(color=mycols, hjust=0.5,vjust = 0.5, size=12, family=fontFamily))

outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withSymb_barplot.", plotType))
ggsave(plot = fract_plot_with_symb, filename = outFile, height=myHeightGG, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))




######################################################################################
# FRACT and FCC AUC - discrete x-axis
######################################################################################

xax_exp <- 0.5
yax_exp <- 0.01


# auc_fract_dt <- get(load("BARPLOT_WITH_FCC_FRACT/all_dt.Rdata"))
# auc_ratio_dt <- get(load("BARPLOT_FCC_AUC_RATIO//all_dt.Rdata"))

auc_fract_dt <- get(load(auc_fract_file))
auc_ratio_dt <- get(load(auc_ratio_file))

auc_fract_ratio_dt <- merge(auc_fract_dt, auc_ratio_dt, by=c("hicds", "exprds"), all.x=TRUE, all.y=TRUE)
stopifnot(!is.na(auc_fract_ratio_dt))

auc_fract_ratio_dt$intervalFCC <- gsub("FCC ", "", auc_fract_ratio_dt$intervalFCC)
# auc_fract_ratio_dt$intervalFCC <- factor(auc_fract_ratio_dt$intervalFCC, levels=rev(fcc_fract_names))
auc_fract_ratio_dt$intervalFCC <- factor(auc_fract_ratio_dt$intervalFCC, levels=rev(gsub("FCC ", "", fcc_fract_names)))

ds_rank_dt <- auc_fract_ratio_dt[,c("hicds", "exprds", "fcc_auc")]
ds_rank_dt$dataset <- paste0(ds_rank_dt$hicds,"\n", ds_rank_dt$exprds)
ds_rank_dt <- unique(ds_rank_dt)
stopifnot(!duplicated(ds_rank_dt$dataset))
ds_rank_dt <- ds_rank_dt[order(ds_rank_dt$fcc_auc, decreasing=TRUE),]
ds_auc_order <- ds_rank_dt$dataset
mycols_scat <- all_cols[all_cmps[paste0(ds_rank_dt$exprds)]]

auc_fract_ratio_dt$dataset <- paste0(auc_fract_ratio_dt$hicds, "\n", auc_fract_ratio_dt$exprds)
auc_fract_ratio_dt$dataset <- factor(auc_fract_ratio_dt$dataset, levels=ds_auc_order)
stopifnot(!is.na(auc_fract_ratio_dt$dataset))

auc_fract_ratio_dt$ds_rank <- as.numeric(auc_fract_ratio_dt$dataset)
xbreaks <- sort(unique(auc_fract_ratio_dt$ds_rank))


myPals <-  eval(parse(text=paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(unique(auc_fract_ratio_dt$intervalFCC)))

scatPlot_tmp <- ggscatter(auc_fract_ratio_dt, 
                      title = paste0("all datasets (n=", length(unique(file.path(auc_fract_ratio_dt$hicds, auc_fract_ratio_dt$exprds))), ")"),
                      x = "ds_rank", 
                      y = "countFCC",
                      color = "intervalFCC",
                      xlab = "Datasets (ranked by decreasing FCC AUC ratio)",
                      ylab = "Ratio of TADs",
                      palette = myPals) +
  labs(color=legTitle)+
  
  
  geom_smooth(aes(color = intervalFCC),method = "lm", linetype=2, se=F)+
  
  theme( # Increase size of axis lines
    
    plot.margin = unit(plotMargin, "lines"), # top, right
    
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14, family=fontFamily),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.title = element_text(face="bold", family=fontFamily)
  )


scatPlot <- scatPlot_tmp +
  
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand=c(yax_exp, yax_exp)) +
  
  scale_x_continuous(labels = rep(labsymbol, length(xbreaks)), breaks= xbreaks, expand=c(xax_exp,xax_exp)) + 
  theme( axis.text.x = element_text(color=mycols_scat, hjust=0.5,vjust = 0.5, size=12, family=fontFamily),
         axis.line.x = element_blank()
         ) +

  # geom_text(data=legDT, aes(label = paste0(labsymbol, " ", legDT$legLabel), x = 60, y =c(0, 0.015, 0.03)),
  #           vjust = 0, hjust=0, size=4,
  #           inherit.aes = FALSE, color = legDT$cmpColor) +
  # 
  
  annotate("text", label = paste0(labsymbol, " ", legDT$legLabel), x = max(xbreaks) + 2, y =c(0, 0.015, 0.03),
           vjust = 0, hjust=0, size=4,
           color = legDT$cmpColor) +
  
  
  annotate("segment", x = 1-xax_exp, y = 0-yax_exp, xend = 58+xax_exp, yend = 0-yax_exp) +
  
  coord_cartesian(clip = 'off', expand=F) +
  
  theme(
    legend.position = "right",
    legend.justification = c("right", "top")
  )

outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_nbrSignifs_auc_scatterplot_xcont.", plotType))
ggsave(plot = scatPlot, filename = outFile, height=myHeightGG, width = myWidthGG*1.8)
cat(paste0("... written: ", outFile, "\n"))

xax_exp <- 0.01


scatPlot_withLabs  <- scatPlot_tmp + 
  
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  
  scale_x_continuous(labels = ds_auc_order, breaks=xbreaks, expand=c(xax_exp,xax_exp))+
  theme(axis.text.x = element_text(color=mycols_scat, hjust=1,vjust = 0.5, size=4, family=fontFamily, angle=90))

outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_nbrSignifs_auc_scatterplot_xcont_withLabs.", plotType))
ggsave(plot = scatPlot_withLabs, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


obs_dt <- get(load(obs_ratio_file))
rd_dt <- get(load(rd_ratio_file))

rdType <- unique(gsub(".+_(.+?)_40kb", "\\1", rd_dt$hicds))

rd_dt$hicds <- gsub(paste0("_", rdType, "_40kb"), "_40kb", rd_dt$hicds)
rd_dt$fcc_auc_rd <- rd_dt$fcc_auc
rd_dt$fcc_auc <- NULL
rd_dt$coexpr_auc_rd <- rd_dt$coexpr_auc
rd_dt$coexpr_auc <- NULL

obs_dt$fcc_auc_obs<- obs_dt$fcc_auc
obs_dt$fcc_auc <- NULL
obs_dt$coexpr_auc_obs <- obs_dt$coexpr_auc
obs_dt$coexpr_auc <- NULL

merged_dt <- merge(rd_dt, obs_dt, by=c("hicds", "exprds"), all=TRUE)
stopifnot(!is.na(merged_dt))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder, paste0("all_ds_fcc_auc_randommidposstrict_vs_obs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x=merged_dt$fcc_auc_obs,
  y=merged_dt$fcc_auc_rd,
  xlab="FCC AUC (OBSERVED)",
  ylab=paste0("FCC AUC (", rdType, ")"),
  main = paste0("OBS vs. ", rdType, " - FCC AUC"),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  pch = 16
)
mtext(side=3, text= paste0("all DS - n=", length(unique(file.path(merged_dt$hicds, merged_dt$exprds)))))
addCorr(x=merged_dt$fcc_auc_obs, y=merged_dt$fcc_auc_rd, bty="n", legPos = "topleft")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("all_ds_coexpr_auc_randommidposstrict_vs_obs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x=merged_dt$coexpr_auc_obs,
  y=merged_dt$coexpr_auc_rd,
  xlab="COEXPR AUC (OBSERVED)",
  ylab=paste0("COEXPR AUC (", rdType, ")"),
  main = paste0("OBS vs. ", rdType, " - COEXPR AUC"),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  pch = 16
)
mtext(side=3, text= paste0("all DS - n=", length(unique(file.path(merged_dt$hicds, merged_dt$exprds)))))
addCorr(x=merged_dt$coexpr_auc_obs, y=merged_dt$coexpr_auc_rd, bty="n", legPos = "topleft")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("all_ds_coexpr_vs_fcc_auc_randommidposstrict.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x=merged_dt$fcc_auc_rd,
  y=merged_dt$coexpr_auc_rd,
  xlab=paste0("FCC AUC (", rdType, ")"),
  ylab=paste0("COEXPR AUC (", rdType, ")"),
  main = paste0("COEXPR vs. FCC AUC - ", rdType),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  pch = 16
)
mtext(side=3, text= paste0("all DS - n=", length(unique(file.path(merged_dt$hicds, merged_dt$exprds)))))
addCorr(x=merged_dt$fcc_auc_rd, y=merged_dt$coexpr_auc_rd, bty="n", legPos = "topleft")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))






outFile <- file.path(outFolder, paste0("all_ds_coexpr_vs_fcc_auc_obs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x=merged_dt$fcc_auc_obs,
  y=merged_dt$coexpr_auc_obs,
  xlab="FCC AUC (OBSERVED)",
  ylab="COEXPR AUC (OBSERVED)",
  main = paste0("COEXPR vs. FCC AUC - OBSERVED"),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  pch = 16
)
mtext(side=3, text= paste0("all DS - n=", length(unique(file.path(merged_dt$hicds, merged_dt$exprds)))))
addCorr(x=merged_dt$fcc_auc_obs, y=merged_dt$coexpr_auc_obs, bty="n", legPos = "topleft")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




