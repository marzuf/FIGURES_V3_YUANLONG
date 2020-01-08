
plotType <- "svg"


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")
require(ggsci)

# Rscript nSignifGenes_nSignifTADs.R

outFolder <- file.path("NSIGNIFGENES_NSIGNIFTADS")
dir.create(outFolder, recursive = TRUE)

dotpch <- 19
# segcol <-  "#BEBEBE19"
segcol <- "grey"
dotCex <- 1.1

geneSignifThresh <- 0.05
tadSignifThresh <- 0.01

inDT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
geneDT <- inDT[,c("hicds", "exprds", "entrezID", "adj.P.Val")]
geneDT <- unique(geneDT)
nSignifGenes_dt <- aggregate(adj.P.Val~hicds + exprds, data = geneDT, FUN=function(x) sum(x<=geneSignifThresh))
colnames(nSignifGenes_dt)[colnames(nSignifGenes_dt) == "adj.P.Val"] <- "nSignifGenes"

tadDT <- inDT[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tadDT <- unique(tadDT)
nSignifTADs_dt <- aggregate(tad_adjCombPval~hicds + exprds, data = tadDT, FUN=function(x) sum(x<=tadSignifThresh))
colnames(nSignifTADs_dt)[colnames(nSignifTADs_dt) == "tad_adjCombPval"] <- "nSignifTADs"

nSignif_dt <- merge(nSignifGenes_dt, nSignifTADs_dt, by=c("hicds", "exprds"), all=TRUE)
stopifnot(!is.na(nSignif_dt))

nSignif_dt <- nSignif_dt[order(nSignif_dt$nSignifTADs, decreasing = TRUE),]
nSignif_dt$dataset <- paste0(nSignif_dt$hicds, "\n", nSignif_dt$exprds)

labcols <- all_cols[all_cmps[nSignif_dt$exprds]]

maxTADs <- max(ceiling(nSignif_dt$nSignifTADs/10)*10)
maxGenes <- max(ceiling(nSignif_dt$nSignifGenes/1000)*1000)

nSignif_dt$nSignifTADs_rescaled <- nSignif_dt$nSignifTADs/maxTADs * maxGenes
nSignif_dt$nSignifGenes_rescaled <- nSignif_dt$nSignifGenes/maxGenes * maxTADs

outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withSymb.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2))
dev.control(displaylist="enable")
par(bty="U", family=fontFamily)
par(mar = c(5,5,2,5))
plot(
  NULL,
  # x = 1:nrow(nSignif_dt),
  # y = nSignif_dt$nSignifTADs, 
  # col = tad_signif_col,
  xlim = c(1, nrow(nSignif_dt)),
  ylim=c(0, maxTADs),
  ylab = "# signif. TADs",
  xlab = "",
  pch=dotpch,
  cex = dotCex,
  main = paste0("# signif. features"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  col.lab = tad_signif_col,
  axes = FALSE
)


segments(
  x0= 1:nrow(nSignif_dt),
  y0 = nSignif_dt$nSignifGenes_rescaled,
  x1= 1:nrow(nSignif_dt),
  y1=nSignif_dt$nSignifTADs,
  col = segcol
)

points(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifTADs, 
  col = tad_signif_col,
  cex = dotCex,
  pch=dotpch
)


legend("topleft",
       legend=c(paste0("gene signif.: p-val <= ", geneSignifThresh),
                paste0("TAD signif.: p-val <= ", tadSignifThresh)),
       bty="n")



mtext(side=3, line=-1, text=paste0("all datasets - n= ", length(unique(nSignif_dt$dataset))))
# mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex =0.5)
axis(2, col = tad_signif_col, col.ticks = tad_signif_col, col.axis=tad_signif_col, at=seq(from=0, to = maxTADs, by=10))
axis(1, labels=F, lwd.ticks = -1)
# axis(1, labels=F, at=1:nrow(nSignif_dt))
par(new = T, family=fontFamily)
plot(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifGenes,
  col = gene_signif_col,
  ylab = NA,
  ylim=c(0, maxGenes),
  xlab = NA,
  pch=dotpch,
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  axes = FALSE
)
axis(side=4, col = gene_signif_col, col.ticks = gene_signif_col, col.axis=gene_signif_col, at = seq(from=0, to=maxGenes, by=1000))
mtext(side = 4, line = 3, '# signif. genes', col=gene_signif_col,  cex=plotCex)



legend("bottom", 
       # legend = paste0(labsymbol, " ", names(all_cols)),
       legend = paste0(names(all_cols)),
       col=all_cols,
       pch=dotpch,
       # lty=c(1,2),
       horiz=TRUE,
       inset=c(0,-0.12),
       cex = plotCex,
       xpd=TRUE, bty="n"
)
signifPlot <- recordPlot()

mtext(side=1, col = labcols, text = rep(labsymbol, nrow(nSignif_dt)), at= 1:nrow(nSignif_dt), las=2, cex = 1.2)


invisible(dev.off())
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withLeg.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2.5))
par(mar = c(5,5,2,5))
replayPlot(signifPlot) 
mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex = 0.6)


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))





