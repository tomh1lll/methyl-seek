#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(ggplot2)

comb_file = read.delim(args[1],h=F)
colnames(comb_file) = c("chr","start","end","min_p","n_probes","z_p","z_sidak_p")


chrs <- factor(c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr23','chrX'))


mergefile2 = comb_file[order(comb_file$chr),]

mergefile2$start = as.numeric(as.character(mergefile2$start))
mergefile2$p = as.numeric(as.character(mergefile2$z_sidak_p))
mergefile2$Log10P = -1*log10(mergefile2$p)
mergefile2$mb = mergefile2$start/1000000
mergefile2$chr <- factor(mergefile2$chr, levels = chrs)

cutoff = -1*(log10(0.05/length(mergefile2$chr)))
ylim <- abs(floor(log10(min(na.omit(mergefile2$p))))) + 2
nCHR <- length(unique(mergefile2$chr))

p1 = ggplot(mergefile2,aes(x=mb,y=Log10P,col=as.factor(chr),size=Log10P)) + geom_point(alpha=0.75) + theme_bw() +
facet_grid(cols = vars(chr),space="free_x",scales="free_x") +
geom_hline(yintercept=cutoff, linetype="dashed") +
labs(x = "BP (Mb)", y = "-log10(p)") +
scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
scale_size_continuous(range = c(0.5,3)) +
scale_x_continuous(breaks=seq(0,300,60)) +
theme(legend.position = "none")

ggsave(args[2], plot = p1,width = 40,height = 10, units = "cm",device = "png")
