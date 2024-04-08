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
  make_option(c("-t", "--type"), type="character", default="pdf",
              help="Output file", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

forward <- opt$input_for
reverse <- opt$input_rev
output <- opt$output
type <- opt$type

forward <- read.csv(forward, header = FALSE, sep = "\t", col.names = c("chrom","start","end","coverage"))
reverse <- read.csv(reverse, header = FALSE, sep = "\t", col.names = c("chrom","start","end","coverage"))
forward$start <- forward$start + 1
forward$end <- forward$end + 1
reverse$start <- reverse$start + 1
reverse$end <- reverse$end + 1

forward.vec <- fill_coverage(forward)
reverse.vec <- fill_coverage(reverse)
reverse.vec <- reverse.vec * -1

df_cov <- data.frame(position=1:length(forward.vec),
                     forward=forward.vec,
                     reverse=reverse.vec)

p <- ggplot(data = df_cov) +
  geom_line(aes(x=position,y=forward)) +
  geom_line(aes(x=position,y=reverse))

ggsave(paste(output,".pdf",sep=''),
       plot = p)


position_percentages <- (df_cov$position / length(df_cov$position))*1000

df_cov_percent <- data.frame(position=1:1000)

df_cov_percent$forward <- unlist(lapply(df_cov_percent$position, function(x) {
  tmp_pos <- max(which(abs(position_percentages - x)==min(abs(position_percentages - x))))
  value <- unlist(df_cov[tmp_pos,"forward"])
  return(value)
}))

df_cov_percent$reverse <- unlist(lapply(df_cov_percent$position, function(x) {
  tmp_pos <- max(which(abs(position_percentages - x)==min(abs(position_percentages - x))))
  value <- df_cov[tmp_pos,"reverse"]
  return(value)
}))

df_cov_percent$position <- df_cov_percent$position / 10

write.csv(df_cov,paste(output,"_absolute_values.csv",sep=""),row.names = FALSE)
write.csv(df_cov_percent,paste(output,"_percent_values.csv",sep=""), row.names = FALSE)
  
