# IGV style


# Rscript conserved_regions_viz.R
# Rscript conserved_regions_viz.R norm_vs_tumor
# Rscript conserved_regions_viz.R subtypes
# Rscript conserved_regions_viz.R wt_vs_mut

# Rscript conserved_regions_viz.R <cmpType>

cat("> START ", "conserved_regions_viz.R", "\n")

startTime <- Sys.time()

plotType <- "svg"

synType <- "inferCARs" #inferCARs or proCARs
 # synType <- "proCARs" #inferCARs or proCARs

source("../FIGURES_V2_YUANLONG/settings.R")

require(ggsci)
tad_col <- pal_d3()(3)[1]
gene_col <- pal_d3()(3)[2]
syn_col <- pal_d3()(3)[3]

consAll_col <- "red"
consAbove_col <- "orange"
consBelow_col <- "black"


outFolder <- file.path("CONSERVED_REGIONS_VIZ")
dir.create(outFolder, recursive = TRUE)


args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  cmpType <- ""
  filePrefix <- ""
  cmpTit <- paste0("all")
} else if(length(args) == 1) {
  cmpType <- args[1]  
  filePrefix <- paste0(cmpType, "_")
  cmpTit <- cmpType
}else {
  stop("---error\n")
}
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3
gene_matching_fuse_threshold <- 0.8


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)


result_dt <- get(load(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/CREATE_FINAL_TABLE/all_result_dt.Rdata")))
nDS <- length(unique(file.path(result_dt$hicds, result_dt$exprds)))

inFolder <- file.path(runFolder, "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2")
outFile <- file.path(inFolder, paste0(filePrefix, "conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(outFile))
conserved_dt <- get(load(outFile))
conserved_dt$conserved_region <- as.character(conserved_dt$conserved_region)

outFile <- file.path(inFolder, paste0(filePrefix, "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(outFile))
conserved_list <- get(load(outFile))
nConserved <- lengths(conserved_list)

stopifnot(length(nConserved) == nrow(conserved_dt))

maxConserved <- names(which.max(nConserved))

max_dt <- conserved_dt[conserved_dt$conserved_region == maxConserved,,drop=F]
stopifnot(nrow(max_dt) == 1)

all_max_regions <- unlist(strsplit(max_dt$corresp_tads, split=","))
stopifnot(length(all_max_regions) == nConserved[maxConserved])

all_max_entrez <- unlist(strsplit(max_dt$intersect_genes_entrez, split=","))
stopifnot(all_max_entrez %in% gff_dt$entrezID)


# all_genes_starts_ends <- sapply(all_max_entrez, function(x) {
#   c(start = gff_dt$start[gff_dt$entrezID == x],
#     end = gff_dt$end[gff_dt$entrezID == x],
#     symbol = gff_dt$symbol[gff_dt$entrezID == x]
#     )
# })

nDScons <- length(all_max_regions)

colConsThresh <- ceiling(nDScons/2)

# retrieve all the genes in the corresponding regions

all_symbs <- as.character(unlist(sapply(all_max_regions, function(ds) {
  reg_symb <- unlist(strsplit(x=result_dt$region_genes[result_dt$hicds == dirname(dirname(ds)) & result_dt$exprds == basename(dirname(ds)) & result_dt$region == basename(ds)], split=","))
  stopifnot(length(reg_symb) > 0 )
  reg_symb
})))
symb_count <- setNames(as.numeric(table(all_symbs)), names(table(all_symbs)))

all_genes_starts_ends <- sapply(unique(all_symbs), function(x) {
  nSymb <- as.numeric(symb_count[paste0(x)])
  c(start = gff_dt$start[gff_dt$symbol == x],
    end = gff_dt$end[gff_dt$symbol == x],
    symbol = x,
    count = nSymb,
    col = ifelse(nSymb == nDScons, consAll_col, ifelse(nSymb >= colConsThresh, consAbove_col, consBelow_col))
  )
})

all_genes_starts_ends <- all_genes_starts_ends[,order(as.numeric(all_genes_starts_ends["start",]))]


all_regions_starts_ends <- sapply(all_max_regions, function(x) {
    g2t_dt <- read.delim(file.path(runFolder, dirname(dirname(x)), "genes2tad", "all_assigned_regions.txt"), header=FALSE, stringsAsFactors=FALSE, col.names=c("chromo", "region", "start", "end"))
    stopifnot(sum(g2t_dt$region == basename(x)) == 1)
    g2t_dt$start[g2t_dt$region == basename(x)]
    c(start=g2t_dt$start[g2t_dt$region == basename(x)], end=g2t_dt$end[g2t_dt$region == basename(x)], chromo = g2t_dt$chromo[g2t_dt$region == basename(x)])
  })

all_exprds <- basename(dirname(colnames(all_regions_starts_ends)))
all_regions_starts_ends <- all_regions_starts_ends[,rev(colnames(all_regions_starts_ends)[order(all_exprds)])]

chromo <- unique(as.character(all_regions_starts_ends["chromo",]))
stopifnot(length(chromo) == 1)


dsSpace <- 0.5
geneSpace <- 0.8

ySynt <- 0.1
syntLwd <- 2

tadOffset <- 50000

textOffset <- 0.05

xStart <- min(as.numeric(as.character(all_regions_starts_ends["start",]))) - tadOffset
xEnd <- max(as.numeric(as.character(all_regions_starts_ends["end",]))) + tadOffset

#### retrieve if there is synteny
synt_dt <- get(load(file.path(runFolder, paste0("SIGNIF_CONSERVED_REGIONS_AND_SYNTENY/", synType, "_overlapDT_bp.Rdata"))))


if(maxConserved %in% synt_dt$refID) {
  syntOffset <- 1
  yOffset <- 0.3 + syntOffset
  
  if(synType == "proCARs") {
    block_dt <- read.delim(file.path(runFolder, "procars_orthology_blocks_processsed.txt"), header=TRUE, stringsAsFactors = FALSE)  
    block_dt$genome[block_dt$genome == "homo_sapiens"] <- "hg19"

  } else if(synType == "inferCARs") {
    block_dt <- read.delim(file.path(runFolder, "inferCARs_data/Orthology.Blocks_processed_hg19.txt"), header=TRUE, stringsAsFactors = FALSE)  
    
  } else {
    stop("-----ok\n")
  }
  synmatch_start <- block_dt$start[
    block_dt$genome == "hg19" &
      block_dt$blockID %in% synt_dt$queryID[synt_dt$refID == maxConserved]]
  synmatch_end <- block_dt$end[
    block_dt$genome == "hg19" &
      block_dt$blockID %in% synt_dt$queryID[synt_dt$refID == maxConserved]]
  
} else {
  yOffset <- 0.3   
  synmatch_start <- NULL
  synmatch_end <- NULL
}


dsPos <- seq(from=0, by=dsSpace, length.out=ncol(all_regions_starts_ends))
genePos <- seq(from = max(dsPos) + dsSpace, by=geneSpace, length.out=ncol(all_genes_starts_ends))
# genePos <- seq(from = max(dsPos) + dsSpace, by=0, length.out=ncol(all_genes_starts_ends))  # all genes on same y position
yStart <- min(c(dsPos, genePos)) - yOffset
yEnd <- max(c(dsPos, genePos)) + yOffset

subTit <- paste0("(# DS conserv. = ", nDScons, ")")

outFile <- file.path(outFolder, paste0(maxConserved, "_viz.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
dev.control(displaylist="enable")
initMar <- par()$mar
par(mar=initMar+c(0,10,0,0))
par(family="Heshley")
par(xpd=TRUE)
plot(NULL,
     main = paste0(maxConserved),
     xlim = c(xStart, xEnd),
     # ylim = c(yStart, yEnd),
     ylim = c(0, yEnd),
     xlab = "",
     # xlab = paste0(gsub("chr", "chromosome ", chromo)),
     ylab = "",
     axes = FALSE,
     cex.main = plotCex
     )
mtext(side = 3, text = subTit)
mtext(side = 1, text= paste0(gsub("chr", "chromosome ", chromo)), font = 2, cex = 1, line=2)
axis(1,
     at = unique(sort(c(as.numeric(as.character(all_regions_starts_ends["start",])), as.numeric(as.character(all_regions_starts_ends["end",]))))), 
     cex = 0.6)
# draw the tads
segments(
  x0 = as.numeric(all_regions_starts_ends["start",]),
  y0 = dsPos,
  x1 = as.numeric(all_regions_starts_ends["end",]),
  y1 = dsPos,
  col=tad_col
)
text(
  x = as.numeric(all_regions_starts_ends["start",]),
  y = dsPos + textOffset,
  # labels = colnames(all_regions_starts_ends),
  labels = dirname(colnames(all_regions_starts_ends)),
  # cex = 0.5,
  cex = 0.7,
  pos=2,
  col = tad_col
)
segments(
  x0 = as.numeric(all_genes_starts_ends["start",]),
  y0 = genePos,
  x1 = as.numeric(all_genes_starts_ends["end",]),
  y1 = genePos,
  col=all_genes_starts_ends["col",]
)
text(
  x = 0.5*(as.numeric(all_genes_starts_ends["start",]) + as.numeric(all_genes_starts_ends["end",])),
  y = genePos - 5*textOffset,
  # y = genePos + 0,
  labels = as.character(all_genes_starts_ends["symbol",]),
  cex = 1,
  pos=3,
  col=as.character(all_genes_starts_ends["col",])
)

segments(x0=c(as.numeric(all_genes_starts_ends["start",]), as.numeric(all_genes_starts_ends["end",])),
         y0=0-yOffset,
         x1=c(as.numeric(all_genes_starts_ends["start",]), as.numeric(all_genes_starts_ends["end",])),
         y1 = rep(genePos,2),
         lty=2, 
         col = as.character(all_genes_starts_ends["col",])
         )


legend(
  "bottom",
  horiz = TRUE,
  legend = c(
    paste0("# cons. = ", nDScons), 
    paste0("# cons. >= ", colConsThresh),
    paste0("# cons. < ", colConsThresh)
    ),
  text.col = c(consAll_col, consAbove_col, consBelow_col), 
  inset = c(-0.0, -0.2),
  bty = "n",
  xpd = TRUE
)

vizplot <- recordPlot()
  

invisible(dev.off())
cat(paste0("... written: ", outFile, "\n"))

if( !is.null(synmatch_start)) {
  outFile <- file.path(outFolder, paste0(maxConserved, "_", synType, "_viz.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
  
  replayPlot(vizplot) 
  
  segments(
    x0 = synmatch_start,
    y0=max(genePos) + syntOffset,
    x1=synmatch_end,
    y1=max(genePos) + syntOffset,
    col = syn_col,
    lwd=syntLwd
  )

  text(
    x = min(synmatch_start),
    y = max(genePos) + syntOffset,
    labels = paste0(synType),
    # cex = 0.5,
    cex = 0.7,
    pos=2,
    col = syn_col
  )

  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
} else{
  cat(paste0("... no syntenic block for conserved region\n"))
}