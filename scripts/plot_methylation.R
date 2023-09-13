#---------------------------------------------------------
# Copyright 2015 Ontario Institute for Cancer Research
# Written by Jared Simpson (jared.simpson@oicr.on.ca)
#---------------------------------------------------------
# obtained from:https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html

library(ggplot2)
library(RColorBrewer)
output <- "bul_vs_f5c.pdf"
input <- "bul_vs_f5c.tsv"
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
    
