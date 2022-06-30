# install.packages('BiocManager') library(BiocManager)
# BiocManager::install('GenomicRanges') BiocManager::install('Biostrings')
# BiocManager::install('Rsamtools') BiocManager::install('GenomicAlignments')
# install.packages('ExomeDepth')

library(ExomeDepth)
library(rlang)

{
  data(genes.hg19)
  data(exons.hg19)

  bedfile <- snakemake@input[["bedfile"]]
  my.bam <- snakemake@input[["bam"]]
  load(snakemake@input[["ref_count"]])

  x = load(snakemake@input[["ref_count"]])
  RefCount.mat = get(x)

  # Probes as bed.frame no header
  bed.hg19 <- read.csv(bedfile, header = F, sep = "\t")

  # Create counts dataframe for BAM
  my.counts <- getBamCounts(bed.frame = bed.hg19, bam.files = my.bam, include.chr = F)

  ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], "data.frame")  #Samma??

  if (names(ExomeCount.dafr[1])=="space"){
  	names(ExomeCount.dafr)[1] <- "chromosome"
  }
  if (names(ExomeCount.dafr[5])=="names"){
    names(ExomeCount.dafr)[5] <- "exon"
  }

  # Remove chr from chromosome column
  ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome),
    pattern = "chr", replacement = "")

  # Prepare the main matrix of read count data
  ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = ".bam$")])  #bara antal counts

  # Create the aggregate reference set for this sample
  my.choice <- select.reference.set(test.counts = ExomeCount.mat[, 1], reference.counts = RefCount.mat,
    bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000, n.bins.reduced = 10000)

  my.reference.selected <- apply(X = RefCount.mat[, my.choice$reference.choice,
    drop = FALSE], MAR = 1, FUN = sum)

  message("Now creating the ExomeDepth object")
  all.exons <- new("ExomeDepth", test = ExomeCount.mat[, 1], reference = my.reference.selected,
    formula = "cbind(test, reference) ~ 1")

# Call CNVs
  all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4, chromosome = ExomeCount.dafr$chromosome,
    start = ExomeCount.dafr$start, end = ExomeCount.dafr$end, name = ExomeCount.dafr$exon)


  # Annotate the ExomeDepth object with Exons
  exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
    IRanges::IRanges(start = exons.hg19$start, end = exons.hg19$end), names = exons.hg19$exon)

  all.exons <- AnnotateExtra(x = all.exons, reference.annotation = exons.hg19.GRanges,
    min.overlap = 1e-04, column.name = "exons")

  # Annotate the ExomeDepth object with Genes
  genes.hg19.GRanges <- GenomicRanges::GRanges(seqnames = genes.hg19$chromosome,
    IRanges::IRanges(start = genes.hg19$start, end = genes.hg19$end), names = genes.hg19$name)
  all.exons <- AnnotateExtra(x = all.exons, reference.annotation = genes.hg19.GRanges,
    min.overlap = 1e-04, column.name = "gene")

  save(all.exons, file = snakemake@output[["exon"]])

  # Prepare output data frame
  output_df = data.frame(all.exons@CNV.calls[order(all.exons@CNV.calls$BF, decreasing = TRUE),])

  # Txt file
  write.table(file = snakemake@output[["result"]], x = output_df, row.names = FALSE, sep = "\t")


  # NexusSV .SV.txt file
  nexus <- c("id", "type")
  nexus_df = output_df[nexus]
  nexus_df$type[nexus_df$type == "duplication"] <- "CN Gain"
  nexus_df$type[nexus_df$type == "deletion"] <- "CN Loss"

  names(nexus_df) <- c("Chromosome Region", "Event")

  write.table(nexus_df, file = snakemake@output[["aggregated_result"]], row.names = FALSE, quote = FALSE, sep = "\t")


  # AED file
  keep <- c("chromosome", "start", "end", "gene", "nexons", "reads.ratio", "type",
    "type")
  aed_df = output_df[keep]
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
