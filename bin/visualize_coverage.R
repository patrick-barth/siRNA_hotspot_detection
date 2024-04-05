#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)

fill_coverage <- function(coverage){
  vec <- rep(NA,max(coverage$end))
  for (i in 1:nrow(coverage)) {
    vec[coverage[i,"start"]:coverage[i,"end"]] <- coverage[i,"coverage"]
  }
  return(vec)
}

option_list <- list(
  make_option(c("-f", "--input_for"), type="character", default=NULL,
              help="Bed file with forward coverage", metavar="character"),
  make_option(c("-r", "--input_rev"), type="character", default=NULL,
              help="Bed file with reverse coverage", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="Output file", metavar="character"),
  make_option(c("-t","--type"), type="character",default="png",
              help="File type of output plot", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

forward <- opt$input_for
reverse <- opt$input_rev

# tmp
forward <- "/home/patrick/tmp/workflow-timo/ID22_41_L29_P8_3_S52_L001_R1_001.for.bed"
reverse <- "/home/patrick/tmp/workflow-timo/ID22_41_L29_P8_3_S52_L001_R1_001.rev.bed"
# tmp end

forward <- read.csv(forward, header = FALSE, sep = "\t", col.names = c("chrom","start","end","coverage"))
reverse <- read.csv(reverse, header = FALSE, sep = "\t", col.names = c("chrom","start","end","coverage"))
forward$start <- forward$start + 1
forward$end <- forward$end + 1
reverse$start <- reverse$start + 1
reverse$end <- reverse$end + 1

forward.vec <- fill_coverage(forward)
reverse.vec <- fill_coverage(reverse)
