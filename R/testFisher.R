#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(dplyr))
df = read.table(args[1], header=F)
colnames(df)[1:4] = c('nU','nD','nQ','nDQ')
suppressMessages(library(sesame))
df2 = with(df, sesame:::testEnrichmentFisherN(nD, nQ, nDQ, nU))
df2 = cbind(df[,colnames(df)[!(colnames(df) %in% c("nD","nQ","nDQ"))],drop=F],df2)
options(scipen=10,digits=5)
write.table(format(df2), file="", quote=FALSE, row.names=FALSE, sep="\t")
## suppressMessages(library(readr))
## cat(readr::format_tsv(df2))

