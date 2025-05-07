#NETCAGE CAGEr analysis
#-----------------------------------------------------------------------
# Use a package
library(dplyr)
library(reshape)
library(gridExtra)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(CAGEr)
library(GenomicFeatures)
library(BSgenome)
library(rtracklayer)
library(ggplot2)
library(ChIPseeker)
library(dplyr)
library(ggforce)
library(ggseqlogo)
library(BSgenome.MorexV3.Gatersleben)

# Load BAM files
inputFiles = list.files(path = "path/to/BAM", recursive = T,  pattern = "\\.bam$", full.names = TRUE)
View(inputFiles)
ce_4NET <- CAGEexp( genomeName     = "BSgenome.MorexV3.Gatersleben"
                    , inputFiles     = inputFiles
                    , inputFilesType = "bam"
                    , sampleLabels   = sub( ".bam", "", basename(inputFiles))
)

# For our data "removeFirstG = TRUE", "correctSystematicG = FALSE" worked best
ce_4NET <- getCTSS(ce_4NET, sequencingQualityThreshold = 10,
                   mappingQualityThreshold = 20, removeFirstG = TRUE,
                   correctSystematicG = FALSE, useMulticore = TRUE, nrCores = 1)

CTSStagCountSE(ce_4NET)
coords <- CTSScoordinatesGR(ce_4NET)
df_coords <- as.data.frame(coords)

genome <- BSgenome.MorexV3.Gatersleben
# Get the annotations as a txdb object
txdb <-makeTxDbFromGFF("path/to/gff/Hv_Morex.pgsb.Jul2020.gff3", format="gff3")
txdb

genes <- genes(txdb)
genes$gene_name <- genes$gene_id

#METHOD TO MAKE annotateCTSS WORKING
setMethod("annotateCTSS", c("CAGEexp", "TxDb"), function (object, ranges){
  g <- genes(ranges)
  g$gene_name <- g$gene_id
  annotateCTSS(ce, g)
  CTSScoordinatesGR(ce)$genes      <- ranges2genes(CTSScoordinatesGR(object), g)
  CTSScoordinatesGR(ce)$annotation <- txdb2annot(CTSScoordinatesGR(object), ranges)
  
  annot <- sapply(CTSStagCountDF(ce)
                  , function(X) tapply(X, CTSScoordinatesGR(object)$annotation, sum))
  colData(ce)[levels(CTSScoordinatesGR(ce)$annotation)] <- DataFrame(t(annot))
  
  validObject(ce)
  object
})

#-----------------------------------------------------------
txdb2annot <- function(ranges, annot) {
  findOverlapsBool <- function(A, B) {
    overlap <- findOverlaps(A, B)
    overlap <- as(overlap, "List")
    any(overlap)
  }
  
  classes <- c("promoter", "exon", "intron", "unknown")
  p <- findOverlapsBool(ranges, trim(suppressWarnings(promoters(annot, 500, 500))))
  e <- findOverlapsBool(ranges, exons(annot))
  t <- findOverlapsBool(ranges, transcripts(annot))
  annot <- sapply( 1:length(ranges), function(i) {
    if      (p[i]) {classes[1]}
    else if (e[i]) {classes[2]}
    else if (t[i]) {classes[3]}
    else           {classes[4]}
  })
  
  annot <- factor(annot, levels = classes)
  Rle(annot)
} 
#' ranges2genes
#' 
#' Assign gene symbol(s) to Genomic Ranges.

#-----------------------------------------------------------------
ce_4NET <- annotateCTSS(ce_4NET, genes)

CTSScoordinatesGR(ce_4NET)$annotation <- txdb2annot(CTSScoordinatesGR(ce_4NET), txdb)
annot <- sapply( CTSStagCountDF(ce_4NET)
                 , function(X) tapply(X, CTSScoordinatesGR(ce_4NET)$annotation, sum))

write.csv(annot, file="output/path/annotationCTSS_4NET.txt")
colData(ce_4NET)[levels(CTSScoordinatesGR(ce_4NET)$annotation)] <- DataFrame(t(annot))
validObject(ce_4NET)
#--------------------------------------------------------
##Merging of replicates and naming of samples
ce_4NET <- mergeSamples(ce_4NET, mergeIndex = c(1,1),
                          mergedSampleLabels = c("NETCAGE_4DAG"))
#Normalization to number of sequenced tags
librarySizes(ce_4NET)
plotReverseCumulatives(ce_4NET, fitInRange = c(5, 10000), onePlot = TRUE)
ce_4NET <- normalizeTagCount(ce_4NET, method = "powerLaw", fitInRange = c(5, 10000), alpha = 1.33, T = 10^6)
ce_4NET[["tagCountMatrix"]]
#CTSS clustering
ce_4NET <- clusterCTSS( object = ce_4NET
                        , threshold = 0.1
                        , thresholdIsTpm = TRUE
                        , nrPassThreshold = 1
                        , method = "distclu"
                        , maxDist = 20
                        , removeSingletons = TRUE
                        , keepSingletonsAbove = 5)
ce_4NET <- cumulativeCTSSdistribution(ce_4NET, clusters = "tagClusters", useMulticore = F)
ce_4NET <- quantilePositions(ce_4NET, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
output4 <- tagClustersGR(ce_4NET, sample = "NETCAGE_4DAG", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

output4_anno <- annotatePeak(output4, tssRegion=c(-500, 100),
                             TxDb=txdb, overlap = "all")
output4_df <- as.data.frame(output4_anno)

write.table(output4_df, file="output/path/netcage_dataset.bed", quote=F, sep="\t", row.names=F, col.names=F)




















