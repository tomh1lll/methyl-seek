#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(rGREAT)
library(ggplot2)

input = read.delim(args[1],h=F)

input2 = input[,c(1,2,3)]

job = submitGreatJob(input2,species = "hg38")

tb = getEnrichmentTables(job,download_by = "tsv")

write.table(tb[[1]],file = args[1])
write.table(tb[[2]],file = args[2])

jpeg(file=args[3],width = 12,height = 6,units = "in",res = 300)
res = plotRegionGeneAssociationGraphs(job)
dev.off()

tb2 = head(tb[[2]],25)

enr_plot = ggplot(tb2,aes(size=GeneFoldEnrich,y=Desc,x=-1*log10(BinomBonfP))) + 
  geom_point() + scale_size_continuous(range = c(0.5,4)) + 
  theme_classic() + scale_x_log10() + labs(x="-Log10(Bonf. Corrected p-value)",size="Fold Enrichment") + 
  theme(text=element_text(size = 14))

ggsave(plot = enr_plot,filename = args[4],width = 12,height = 12,units = "in",dpi = "retina")
