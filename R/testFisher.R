#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(dplyr))
options(scipen=99)
options(digits=3)
df = read.table(args[1], header=F)
colnames(df)[1:4] = c('nu','nf','nq','nfq')
df = df %>% dplyr::mutate(
    nfmq = nf - nfq,
    nqmf = nq - nfq,
    numfq = nu - nf - nq + nfq,
    ## https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
    p_val_enr = sprintf("%1.3g", phyper(nfq-1, nf, nu-nf, nq, lower.tail=FALSE)),
    p_val_dep = sprintf("%1.3g", phyper(nfq, nf, nu-nf, nq, lower.tail=TRUE)),
    odds_ratio=(nfq/nfmq)/(nqmf/numfq)) %>% dplyr::arrange(as.numeric(p_val_enr))
write.table(format(df), file="", quote=FALSE, row.names=FALSE, sep="\t")

