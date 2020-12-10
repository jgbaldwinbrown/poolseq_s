#!/usr/bin/env Rscript

suppressMessages(library(poolSeq))
# options(error = recover)

get_info <- function(infopath) {
    info_unstructured = as.data.frame(fread(infopath, sep="\t", header=TRUE))
    info = vector(mode = "list", length = 8)
    names(info) = c("chrom", "pos", "gen", "repl", "gen_levels", "repl_levels", "chrom_levels", "pos_levels")
    info$gen = info_unstructured$gen
    info$repl = info_unstructured$repl
    info$gen_levels = as.numeric(sort(levels(factor(info_unstructured$gen))))
    info$repl_levels = as.numeric(sort(levels(factor(info_unstructured$repl))))
    return(info)
}

update_info <- function(info, sync) {
    info$chrom = sync@alleles$chr
    info$pos = sync@alleles$pos
    info$chrom_levels = sort(levels(factor(info$chrom)))
    info$pos_levels = as.numeric(sort(levels(factor(info$pos))))
    return(info)
}

#        myTraj_repltemp = af.traj(mySync, info$chrom, info$pos, repl)
estimateSH_individual_loci <- function(sync, Ne, gen_levels) {
    sync
}

main = function() {
    
    args = commandArgs(trailingOnly = TRUE)
    syncpath = args[1]
    infopath = args[2]
    info = get_info(infopath)
    mySync <- read.sync(file=syncpath, gen=info$gen, repl=info$repl)
    info = update_info(info, mySync)
    print("mysync")
    print(str(mySync))
    print("info")
    print(str(info))
    myTraj = af.traj(mySync, info$chrom, info$pos, info$repl)
    print("mytraj")
    print(myTraj)
    est_nes = rep(NA, length(info$repl_levels))
    for (repl in info$repl_levels) {
        myTraj_repltemp = af.traj(mySync, info$chrom, info$pos, repl)
        myCov_repltemp = coverage(mySync, info$chrom, info$pos, repl=repl, gen=info$gen_levels)
        traj_gen1_name = paste("F", as.character(info$gen_levels[1]), sep="")
        traj_gen2_name = paste("F", as.character(info$gen_levels[length(info$gen_levels)]), sep="")
        cov_gen1_name = paste("F", as.character(info$gen_levels[1]), ".R", as.character(repl), ".cov", sep="")
        cov_gen2_name = paste("F", as.character(info$gen_levels[length(info$gen_levels)]), ".R", as.character(repl), ".cov", sep="")
        est_nes[repl] = estimateNe(
            p0=myTraj_repltemp[,traj_gen1_name], 
            pt=myTraj_repltemp[,traj_gen2_name], 
            cov0=myCov_repltemp[,cov_gen1_name], 
            covt=myCov_repltemp[,cov_gen2_name], 
            t=info$gen_levels[length(info$gen_levels)] - info$gen_levels[1]
        )
        # est_nes[repl] = estimateNe(
        #     p0=mytraj1[,traj_gen1_name], 
        #     pt=mytraj1[,"F10"], 
        #     cov0=mycov1[,"F0.R1.cov"], 
        #     covt=mycov1[,"F10.R1.cov"], 
        #     t=10
        # )
        # note: add options when ready: Ncensus=1000, poolSize=c(300, 300)
    }
    mean_ne = mean(est_nes)
    # est_ne = estimateNe(p0=mytraj1[,"F0"], pt=mytraj1[,"F10"], cov0=mycov1[,"F0.R1.cov"], covt=mycov1[,"F10.R1.cov"], t=10, Ncensus=1000, poolSize=c(300, 300))
    
    print("est_nes")
    print(est_nes)
    print(str(est_nes))
    
    est_all_p <- estimateSH_windows(myTraj, mean_ne, info$gen_levels)
    est_p <- estimateSH(myTraj, Ne=mean_ne, t=info$gen_levels, h=0.5, simulate.p.value=TRUE)
    print("est_p")
    print(est_p)
    print(str(est_p))
    # print(confint(est))
}

main()
