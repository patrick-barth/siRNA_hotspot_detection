#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(reshape2)

fill_coverage <- function(coverage){
  vec <- rep(NA,max(coverage$end))
  for (i in 1:nrow(coverage)) {
    vec[coverage[i,"start"]:coverage[i,"end"]] <- coverage[i,"coverage"]
  }
  return(vec)
}

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


input <- read.csv(input, header = TRUE, row.names = 1,sep = "\t")

input$row <- seq_len(nrow(input))
dat <- melt(input,id.vars = 'row')

p <- ggplot(dat,aes(x=row,y=value,fill=variable)) +
      geom_bar(stat='identity')

ggsave(paste(output,".pdf",sep=''),
       plot = p)