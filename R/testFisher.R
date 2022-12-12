#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(dplyr))
options(scipen=99)
options(digits=3)
df = read.table(args[1], header=F)
colnames(df)[1:4] = c('nU','nD','nQ','nDQ')
suppressMessages(library(sesame))
df2 = with(df, sesame:::testEnrichmentFisherN(nD, nQ, nDQ, nU))
df2 = cbind(df[,colnames(df)[!(colnames(df) %in% c("nD","nQ","nDQ"))],drop=F],df2)
write.table(format(df2), file="", quote=FALSE, row.names=FALSE, sep="\t")

