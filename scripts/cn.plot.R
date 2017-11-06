#!/usr/bin/env Rscript
#==============================================================================
# title           : cn.plot.R
# description     : Plot allele-specific copy number states from BED files.
# author          : Harald Detering (harald.detering@gmail.com)
# date            : 2017-10-27
# version         : 0.1    
# usage           : Rscript --vanilla cn.plot.R /path/to/bed [filename_pattern]
# notes           :
# R_version       : 3.4.2
#==============================================================================

suppressMessages(require(readr))    # for read_tsv()
suppressMessages(require(dplyr))    # for mutate()
suppressMessages(require(tidyr))    # for unnest()
suppressMessages(require(purrr))    # for map(), reduce()
suppressMessages(require(stringr))  # for gsub()
suppressMessages(require(reshape2)) # for melt()
suppressMessages(require(ggplot2))

data.dir = "."
file.pattern = "*.cn.bed"

# read command line params
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  data.dir <- args[1]
  if (!dir.exists(data.dir)) {
    stop(sprintf("Directory does not exist: %s", data.dir), call.=FALSE)
  }
}
if (length(args) > 1) {
  file.pattern <- args[2]
}

#data.dir = "~/sequenza/c4.r4.m6000.60Mb.100x"

# get input files
bed.files <- dir(data.dir, pattern=file.pattern)
if (length(bed.files) == 0) {
  stop(sprintf("No input files found for pattern: '%s'", file.pattern), call.=FALSE)
}

# get copy number states from files
df.cn <- data_frame(sample = gsub("\\.cn\\.bed", "", bed.files)) %>%
  mutate(file_data = map(bed.files, 
                         ~ read_tsv(file.path(data.dir, .), 
                                    col_names=c("chr","start","end","A","B"))))
df.cn <- unnest(df.cn)

# plot allele-specific copy number
df <- melt(df.cn, id.vars=c("sample","chr","start","end"), variable.name="allele", value.name="copies")
df[df$allele=="A",]$copies <- df[df$allele=="A",]$copies - 0.1
df[df$allele=="B",]$copies <- df[df$allele=="B",]$copies + 0.1

pdf(file.path(data.dir, "cn.true.pdf"))
ggplot(df) + geom_segment(aes(x=start, xend=end, y=copies, yend=copies, color=allele)) + facet_grid(sample ~ chr) +
  theme_bw()
dev.off()

