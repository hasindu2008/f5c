#! /usr/bin/env Rscript

#---------------------------------------------------------
# Copyright 2015 Ontario Institute for Cancer Research
# Written by Jared Simpson (jared.simpson@oicr.on.ca)
#---------------------------------------------------------
# obtained from:https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html
# Usage: Rscript --vanilla plot_methylation.R -i bi_vs_f5c.tsv -o bi_vs_f5c.pdf

#install.packages("optparse")
#install.packages("ggplot2")

library(ggplot2)
library(RColorBrewer)
library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
              help="input methylation comparison tsv file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.pdf",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least -i input.tsv must be provided.n", call.=FALSE)
}

output <- opt$out
input <- opt$input
pdf(file = output, width=10, height=8)
data <- read.table(input, header=T)

# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(data$frequency_1, data$frequency_2)
title <- sprintf("N = %d r = %.10f", nrow(data), c)
ggplot(data, aes(frequency_1, frequency_2)) +
  geom_bin2d(bins=25) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Bisulfite Methylation Frequency") +
  ylab("f5c Methylation Frequency") +
  theme_bw(base_size=20) +
  ggtitle(title)

dev.off()

