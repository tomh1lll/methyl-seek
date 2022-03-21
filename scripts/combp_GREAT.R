#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(rGREAT)

input = read.delim(args[1],h=F)

input2 = input[,c(1,2,3)]

job = submitGreatJob(input2,species = "hg38")

tb = getEnrichmentTables(job,download_by = "tsv")

write.table(tb[[1]],file = args[3])
write.table(tb[[2]],file = args[4])

jpeg(file=args[2],width = 12,height = 6,units = "in",res = 300)
res = plotRegionGeneAssociationGraphs(job)
dev.off()
