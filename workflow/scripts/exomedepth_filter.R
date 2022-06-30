# install.packages('BiocManager') library(BiocManager)
# BiocManager::install('GenomicRanges') BiocManager::install('Biostrings')
# BiocManager::install('Rsamtools') BiocManager::install('GenomicAlignments')
# install.packages('ExomeDepth')

library(ExomeDepth)
library(rlang)

{
  #VARIANT FILES
  Conifer <- read.csv("snakemake@input[["conifer"]]", header=TRUE, sep = "\t")
  data(Conrad.hg19)
  ED <- read.csv("snakemake@input[["ED_common"]]", header=TRUE, sep = "\t")
  load(snakemake@input[["exon"]])
  WW25 <- read.csv("snakemake@input[["WW25"]]", header=TRUE, sep = "\t")


  #CONRAD
  all.exons <- AnnotateExtra(x = all.exons,
                  reference.annotation = Conrad.hg19.common.CNVs,
                  min.overlap = 0.8,
                  column.name = 'Conrad')

  #WW25
  WW25.GRanges <- GenomicRanges::GRanges(seqnames = WW25$chromosome,
                  IRanges::IRanges(start=WW25$start,end=WW25$end),
                  names = WW25$name)

  all.exons <- AnnotateExtra(x = all.exons,
                  reference.annotation = WW25.GRanges,
                  min.overlap = 0.8,
                  column.name = 'WW25')

  #Conifer
  Conifer.GRanges <- GenomicRanges::GRanges(seqnames = Conifer$chromosome,
                  IRanges::IRanges(start=Conifer$start,end=Conifer$end),
                  names = Conifer$name)

  all.exons <- AnnotateExtra(x = all.exons,
                  reference.annotation = Conifer.GRanges,
                  min.overlap = 0.8,
                  column.name = 'Conifer')

  Conifer.GRanges <- GenomicRanges::GRanges(seqnames = Conifer$chromosome,
                  IRanges::IRanges(start=Conifer$start,end=Conifer$end),
                  names = Conifer$sample)

  all.exons <- AnnotateExtra(x = all.exons,
                  reference.annotation = Conifer.GRanges,
                  min.overlap = 0.8,
                  column.name = 'Sample_ConIFER')

  #ExomeDepth
  ED.GRanges <- GenomicRanges::GRanges(seqnames = ED$chromosome,
                IRanges::IRanges(start=ED$start,end=ED$end), names = ED$id)

  all.exons <- AnnotateExtra(x = all.exons, reference.annotation = ED.GRanges,
                min.overlap = 0.8, column.name = 'EDcommon')
​
  ED.GRanges <- GenomicRanges::GRanges(seqnames = ED$chromosome,
              IRanges::IRanges(start=ED$start,end=ED$end), names = ED$CNV_type)

  all.exons <- AnnotateExtra(x = all.exons,
              reference.annotation = ED.GRanges, min.overlap = 0.8,
              column.name = 'CNV_type_ED')


  #NexusSV file
  df = data.frame(all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),])
  filt_df <- df[is.na(df$EDcommon),]
  nexus <- c("id", "type")
  nexus_df = filt_df[nexus]
  nexus_df$type[nexus_df$type == "duplication"] <- "CN Gain"
  nexus_df$type[nexus_df$type == "deletion"] <- "CN Loss"

  names(nexus_df) <- c("Chromosome Region", "Event")

  write.table(nexus_df, file = snakemake@output[["aggregated_result"]], row.names = FALSE, quote = FALSE, sep = "\t")


  # AED file
  keep <- c("chromosome", "start", "end", "gene", "nexons", "reads.ratio", "type",
    "type")
  aed_df = filt_df[keep]
  aed_df$chromosome <- sub("^", "chr", aed_df$chromosome)
  aed_df$type[aed_df$type == "duplication"] <- "copynumber/gain"
  aed_df$type[aed_df$type == "deletion"] <- "copynumber/loss"
  aed_df$type.1[aed_df$type.1 == "duplication"] <- "rgb(0,0,255)"
  aed_df$type.1[aed_df$type.1 == "deletion"] <- "rgb(255,0,0)"

  aed_df$new <- NA  # blank column
  NewNames <- c("bio:sequence(aed:String)", "bio:start(aed:Integer)", "bio:end(aed:Integer)",
    "aed:name(aed:String)", "bio:markerCount(aed:Integer)", "bio:state(aed:Rational)",
    "aed:category(aed:String)", "style:color(aed:Color)", "aed:value(aed:String)")
  names(aed_df) <- NewNames
  aed_df <- aed_df[, c("bio:sequence(aed:String)", "bio:start(aed:Integer)", "bio:end(aed:Integer)",
    "aed:name(aed:String)", "aed:value(aed:String)", "bio:markerCount(aed:Integer)",
    "bio:state(aed:Rational)", "aed:category(aed:String)", "style:color(aed:Color)")]

  header2 <- data.frame("", "", "", "affx:ucscGenomeVersion(aed:String)", "hg19",
    "", "", "", "", stringsAsFactors = FALSE)
  names(header2) <- c("bio:sequence(aed:String)", "bio:start(aed:Integer)", "bio:end(aed:Integer)",
    "aed:name(aed:String)", "aed:value(aed:String)", "bio:markerCount(aed:Integer)",
    "bio:state(aed:Rational)", "aed:category(aed:String)", "style:color(aed:Color)")
  aed_df <- rbind(header2, aed_df)
  header1 <- data.frame("", "", "", "namespace:affx(aed:URI)", "http://affymetrix.com/ontology/",
    "", "", "", "", stringsAsFactors = FALSE)
  names(header1) <- c("bio:sequence(aed:String)", "bio:start(aed:Integer)", "bio:end(aed:Integer)",
    "aed:name(aed:String)", "aed:value(aed:String)", "bio:markerCount(aed:Integer)",
    "bio:state(aed:Rational)", "aed:category(aed:String)", "style:color(aed:Color)")
  aed_df <- rbind(header1, aed_df)

  write.table(aed_df, file = snakemake@output[["aed"]], row.names = FALSE, quote = FALSE, sep = "\t")

}
