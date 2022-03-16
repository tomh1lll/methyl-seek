#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(rGREAT)

input = read.delim(args[1],h=F)

input2 = input[,c(1,2,3)]

job = submitGreatJob(input2)

tb = getEnrichmentTables(job)

write.table(tv[[1]],file = args[3])
write.table(tv[[2]],file = args[4])

jpeg(file=args[2],width = 12,height = 6,units = "in",res = 300)
res = plotRegionGeneAssociationGraphs(job)
dev.off()