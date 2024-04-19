#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="File containing percentages", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="Output file", metavar="character"),
  make_option(c("-t", "--type"), type="character", default="pdf",
              help="Output file", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input <- opt$input
output <- opt$output
type <- opt$type

input <- read.csv(input, header = TRUE,sep = "\t")

p <- ggplot(input[1:11,],aes(x=length,y=count)) +
  geom_bar(stat='identity')

ggsave(paste(output,".pdf",sep=''),
       plot = p)