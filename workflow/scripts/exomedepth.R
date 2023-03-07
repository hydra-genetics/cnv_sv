# install.packages('BiocManager') library(BiocManager)
# BiocManager::install('GenomicRanges') BiocManager::install('Biostrings')
# BiocManager::install('Rsamtools') BiocManager::install('GenomicAlignments')
# install.packages('ExomeDepth')

library(ExomeDepth)

bedfile <- snakemake@input[["bedfile"]]
my.bam <- snakemake@input[["bam"]]
load(snakemake@params[["ref_count"]]) # load refcount_df

message(paste('Exomedepth reference:', snakemake@params[["ref_count"]]))

genome_version <- snakemake@config[["exomedepth_call"]][["genome_version"]]
if (genome_version == "hg38") {

  exons <-  read.csv(snakemake@input[["exons"]],
                     header = TRUE, sep = "\t")
  genes <-  read.csv(snakemake@input[["genes"]],
                     header = TRUE, sep = "\t")

} else  if (genome_version == "hg19") {

  data(genes.hg19) # use genes and exons info packaged with ExomeDepth
  data(exons.hg19)

  exons <- exons.hg19 # rename so it works for code below
  gene <- genes.hg19
  rm(exons.hg19, genes.hg19)

} else {

  stop("Please specify genome_version in the config file 
      as 'hg19' or 'hg38' for the exomedepth_call rule")

}

# Probes as bed.frame no header
targets_bed <- read.csv(bedfile, header = FALSE, sep = "\t")

# Create counts dataframe from the BAM file
ExomeCount.dafr  <- getBamCounts(bed.frame = targets_bed, bam.files = my.bam,
                                 include.chr = FALSE)
save.image('test.Rdata')
## Check that the ExomeCount.dafr and reference df are ordered exactly the same
## by comparing chromosome, start and end. Coordinates are 1-based in both
if (!all.equal(ExomeCount.dafr[,1:3], refcount_df[, 1:3])) {
  stop("The reference used has a differenct sort order to the test sample")
}

# covert the loaded reference count df to a matrix
RefCount.mat <- as.matrix(
  refcount_df[, grep(names(refcount_df), pattern = ".bam$")])

# Remove chr from chromosome column
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome),
  pattern = "chr", replacement = "")

# Prepare the main matrix of read count data
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr),
 pattern = ".bam$")])

# Create the aggregate reference set for this sample
my.choice <- select.reference.set(
  test.counts = ExomeCount.mat[, 1],
  reference.counts = RefCount.mat,
  bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start) / 1000,
  n.bins.reduced = 10000)

my.reference.selected <- apply(
  X = RefCount.mat[, my.choice$reference.choice, drop = FALSE],
  MAR = 1, FUN = sum)

message("Now creating the ExomeDepth object")
all.exons <- new("ExomeDepth", test = ExomeCount.mat[, 1],
                 reference = my.reference.selected,
                 formula = "cbind(test, reference) ~ 1")

message("Now calling CNVs")
all.exons <- CallCNVs(x = all.exons, 
                      transition.probability = 10^-4,
                      chromosome = ExomeCount.dafr$chromosome,
                      start = ExomeCount.dafr$start, end = ExomeCount.dafr$end,
                      name = ExomeCount.dafr$exon)

if (length(all.exons@CNV.calls) > 0) {
  # Annotate the ExomeDepth object with Exons
  message("Annotating with exons")
  exons.GRanges <- GenomicRanges::GRanges(seqnames = exons$chromosome,
    IRanges::IRanges(start = exons$start, end = exons$end), names = exons$name)

  all.exons <- AnnotateExtra(x = all.exons, 
                             reference.annotation = exons.GRanges,
                             min.overlap = 1e-04, column.name = "exons")

  # Annotate the ExomeDepth object with Genes
  message('Annotating with genes')
  genes.GRanges <- GenomicRanges::GRanges(
    seqnames = genes$chromosome,
    IRanges::IRanges(start = genes$start, end = genes$end),
    names = genes$name)

  all.exons <- AnnotateExtra(
    x = all.exons,
    reference.annotation = genes.GRanges,
    min.overlap = 1e-04, 
    column.name = "gene")

  # save the cnv calls object as Rdata file
  save(all.exons, file = snakemake@output[["exon"]])

  # Prepare output data frame
  output_df <- data.frame(all.exons@CNV.calls[order(all.exons@CNV.calls$BF,
                                                   decreasing = TRUE),])

  # Txt file
  write.table(file = snakemake@output[["txt"]], x = output_df,
              row.names = FALSE, sep = "\t")

} else {
  message("No result found")
  writeLines("", snakemake@output[["txt"]])
  writeLines("", snakemake@output[["exon"]])
}