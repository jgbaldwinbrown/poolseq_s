#!/usr/bin/env Rscript

suppressMessages(library(Nest))
suppressMessages(library(data.table))

################################################################################
# functions:

ne_pair <- function(data, poptable, popgroup, gen0, gen1) {
    pop0 = poptable[poptable$group == popgroup & poptable$time_point==gen0,]
    pop1 = poptable[poptable$group == popgroup & poptable$time_point==gen1,]
    p0 = data[,pop0$afname[1]]
    p1 = data[,pop1$afname[1]]
    cov0 = data[,pop0$covname[1]]
    cov1 = data[,pop1$covname[1]]
    gendiff = pop1$time_point[1] - pop0$time_point[1]
    return(c(popgroup, gen0, gen1, estimateNe(p0, p1, cov0, cov1, gendiff, Ncensus=50, method="P.planII")))
}

ne_whole_group <- function(data, poptable, popgroup) {
    minitable <- poptable[poptable$group == popgroup,]
    gens = sort(minitable$time_point)
    ngens = length(gens)
    nout = (ngens * (ngens - 1)) / 2
    out = data.frame(group = rep(NA, nout), gen0 = rep(NA, nout), gen1 = rep(NA, nout), ne = rep(NA, nout))
    outi = 1
    
    for (i in 1:(length(gens)-1)) {
        gen0 = gens[i]
        for (j in (i+1):length(gens)) {
            gen1 = gens[j]
            out[outi,] = ne_pair(data, poptable, popgroup, gen0, gen1)
            outi = outi+1
        }
    }
    return(out)
}

ne_all_groups <- function(data, poptable) {
    groups = sort(levels(factor(poptable$group)))
    ngroups = length(groups)
    out = ne_whole_group(data, poptable, groups[1])
    if (ngroups >= 2) {
        for (i in 2:ngroups) {
            out_temp = ne_whole_group(data, poptable, groups[i])
            out = as.data.frame(rbind(out, out_temp))
        }
    }
    return(out)
}

################################################################################

main <- function() {
    args = commandArgs(trailingOnly=TRUE)
    datapath = args[1]
    poptablepath = args[2]
    data = as.data.frame(fread(datapath, header=TRUE))
    poptable = as.data.frame(fread(poptablepath, header=TRUE))
    out = ne_all_groups(data, poptable)
    write.table(out, "", sep="\t", quote=FALSE)
}

main()
