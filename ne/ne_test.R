#!/usr/bin/env Rscript

library(Nest)
library(data.table)

data = as.data.frame(fread("temp.txt", header=TRUE))
str(data)
p0 = data$NA00001_af
pt = data$NA00002_af
cov0 = data$NA00001_count
covt = data$NA00002_count

# $ NA00001_af   : num  0.0995 0.0995 0.0994 0.0995 0.0995 ...
# $ NA00001_count: int  1005 1085 805 1005 1085 805 1005 1085 805
# $ NA00002_af   : num  0.00991 0.00991 0.00991 0.00991 0.00991 ...
# $ NA00002_count: int  5045 5045 5045 5045 5045 5045 5045 5045 5045
# $ NA00003_af   : num  0.000989 0.000989 0.000989 0.000989 0.000989 ...
# $ NA00003_count: int  10108 10108 10108 10108 10108 10108 10108 10108 10108
# $ NA00004_af   : num  0.00992 0.00992 0.00992 0.00992 0.00992 ...
# $ NA00004_count: int  5545 5545 5545 5545 5545 5545 5545 5545 5545

nes = estimateNe(p0, pt, cov0, covt, 36, Ncensus=50, method="P.planII")
str(nes)
