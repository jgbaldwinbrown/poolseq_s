#!/usr/bin/env Rscript

suppressMessages(library(poolSeq))
suppressMessages(library(parallel))
# options(error = recover)

global_chunksize = 1000

get_info <- function(infopath) {
    info_unstructured = as.data.frame(fread(infopath, sep="\t", header=TRUE))
    # print(info_unstructured)
    info = vector(mode = "list", length = 10)
    names(info) = c("chrom", "pos", "gen", "repl", "pool_size", "gen_levels", "repl_levels", "pool_size_levels", "chrom_levels", "pos_levels")
    info$gen = info_unstructured$gen
    info$repl = info_unstructured$repl
    info$pool_size = info_unstructured$pool_size
    info$gen_levels = sort(as.numeric(levels(factor(info_unstructured$gen))))
    info$repl_levels = sort(as.numeric(levels(factor(info_unstructured$repl))))
    info$pool_size_levels = sort(as.numeric(levels(factor(info_unstructured$pool_size))))
    # print(info)
    return(info)
}

update_info <- function(info, sync) {
    info$chrom = sync@alleles$chr
    info$pos = sync@alleles$pos
    info$chrom_levels = sort(levels(factor(info$chrom)))
    info$pos_levels = as.numeric(sort(levels(factor(info$pos))))
    return(info)
}

estimateSH_one_locus = function(sync, Ne, info, chrom, pos) {
    # print("sync: ")
    # print(sync)
    # print("chrom: ")
    # print(chrom)
    # print("pos: ")
    # print(pos)
    # print("info$repl: ")
    # print(info$repl)
    # traj = af.traj(sync, chrom, pos, repl=info$repl) # TESTING
    traj = af.traj(sync, chrom, pos, repl=info$repl_levels)
    # print("traj: ")
    # print(traj)
    # print("Ne: ")
    # print("info$gen_levels: ")
    est_p <- estimateSH(traj, Ne=round(Ne), t=info$gen_levels, h=0.5, simulate.p.value=TRUE)
    # print(est_p)
    # print(str(est_p))
    # print("est_p: ")
    # print(est_p)
    return(est_p)
}

#        myTraj_repltemp = af.traj(mySync, info$chrom, info$pos, repl)
estimateSH_individual_loci <- function(sync, Ne, info) {
    chrom_pos = cbind(info$chrom, info$pos)
    chrom_pos_list = split(chrom_pos, seq(nrow(chrom_pos)))
    # print("chrom_pos_list: ")
    # print(chrom_pos_list)
    # print("sync: ")
    # print(sync)
    # print("Ne: ")
    # print(Ne)
    # print("info: ")
    # print(info)
    s_vals = unname(mclapply(chrom_pos_list, function(x) {estimateSH_one_locus(sync, Ne, info, x[1], x[2])}))
    # print(s_vals)
    return(list(chrom_pos, s_vals))
}

combine_est_sh_outputs <- function(full_output_list) {
    # print("full_output_list:")
    # print(full_output_list)
    chrom_pos_list = mclapply(full_output_list, function(x) x[[1]])
    s_vals_list = unlist(mclapply(full_output_list, function(x) x[[2]]), recursive=FALSE)
    chrom_pos = do.call("rbind", chrom_pos_list)
    # print(full_output_list)
    # print(chrom_pos_list)
    # print("s_vals_list:")
    # print(s_vals_list)
    # print("s_vals_list[[1]]:")
    # print(s_vals_list[[1]])
    # print(chrom_pos)
    # print(list(chrom_pos, s_vals_list))
    return(list(chrom_pos, s_vals_list))
}

estimateSH_individual_loci_savewrapper <- function(sync, Ne, info, outpath) {
    iteration_series = seq(1,length(info$chrom),global_chunksize)
    full_output_list = vector(mode = "list", length = length(iteration_series))
    j = 1
    for (i in iteration_series) {
        temppath_est_sh = paste(outpath, "_tempdir/", outpath, "_est_sh_", as.character(i), ".RData", sep="")
        temppath_est_sh_done = paste(outpath, "_tempdir/", outpath, "_est_sh_", as.character(i), ".RData.done", sep="")
        if (! file.exists(temppath_est_sh_done)) {
            mini_info = info
            mini_info$chrom = info$chrom[i:min((i+global_chunksize-1), length(info$chrom))]
            mini_info$pos = info$pos[i:min((i+global_chunksize-1), length(info$pos))]
            temp = estimateSH_individual_loci(sync, Ne, mini_info)
            full_output_list[[j]] = temp
            saveRDS(temp, file = temppath_est_sh)
            # print("temp:")
            # print(str(temp))
            # print("temp[[1]]:")
            # print(str(temp[[1]]))
            # print("temp[[2]]:")
            # print(str(temp[[2]]))
            file.create(temppath_est_sh_done)
        } else {
            full_output_list[[j]] = readRDS(temppath_est_sh)
        }
        j = j + 1;
    }
    # print(full_output_list)
    return(combine_est_sh_outputs(full_output_list))
}

gets_from_ests <- function(ests) {
    s_vals = mclapply(ests, function(x) {x$s})
    return(unlist(s_vals))
}

getp_from_ests <- function(ests) {
    p_vals = mclapply(ests, function(x) {x$p.value})
    return(unlist(p_vals))
}

est_full <- function(sync, Ne, info) {
    est_list <- estimateSH_individual_loci(sync, Ne, info)
    # print(est_list)
    chrom_pos = est_list[[1]]
    est_all_p = est_list[[2]]
    s_vals = gets_from_ests(est_all_p)
    p_vals = getp_from_ests(est_all_p)
    out = as.data.frame(cbind(chrom_pos, s_vals, p_vals))
    colnames(out) = c("chrom", "pos", "s", "p.value")
    return(out)
}

est_full_save <- function(sync, Ne, info, outpath) {
    est_list <- estimateSH_individual_loci_savewrapper(sync, Ne, info, outpath)
    # print(est_list)
    chrom_pos = est_list[[1]]
    est_all_p = est_list[[2]]
    s_vals = gets_from_ests(est_all_p)
    p_vals = getp_from_ests(est_all_p)
    # print(s_vals)
    # print(p_vals)
    out = as.data.frame(cbind(chrom_pos, s_vals, p_vals))
    # print(out)
    colnames(out) = c("chrom", "pos", "s", "p.value")
    return(out)
}

main = function() {
    
    args = commandArgs(trailingOnly = TRUE)
    syncpath = args[1]
    infopath = args[2]
    outpath = args[3]
    info = get_info(infopath)
    mySync <- read.sync(file=syncpath, gen=info$gen, repl=info$repl)
    info = update_info(info, mySync)
    # myTraj = af.traj(mySync, info$chrom, info$pos, info$repl)
    
    outdir = paste(outpath, "_tempdir/", sep="")
    if (! dir.exists(outdir)) {
        dir.create(outdir)
    }
    
    est_nes = rep(NA, length(info$repl_levels))
    for (repl in info$repl_levels) {
        myTraj_repltemp = af.traj(mySync, info$chrom, info$pos, repl)
        myCov_repltemp = coverage(mySync, info$chrom, info$pos, repl=repl, gen=info$gen_levels)
        # print(myTraj_repltemp)
        # print(myCov_repltemp)
        gen1 = info$gen_levels[1]
        gen2 = info$gen_levels[length(info$gen_levels)]
        # print(info$gen == gen1)
        # print(info$gen == gen2)
        # print(info$repl == repl)
        # print(info$pool_size)
        pool1 = info$pool_size[info$gen == gen1 & info$repl == repl]
        pool2 = info$pool_size[info$gen == gen2 & info$repl == repl]
        traj_gen1_name = paste("F", as.character(gen1), sep="")
        traj_gen2_name = paste("F", as.character(gen2), sep="")
        cov_gen1_name = paste("F", as.character(gen1), ".R", as.character(repl), ".cov", sep="")
        cov_gen2_name = paste("F", as.character(gen2), ".R", as.character(repl), ".cov", sep="")
        repl_ne_outpath = paste(outdir, outpath, "_", cov_gen1_name, "_ne.txt", sep="")
        repl_ne_outpath_done = paste(outdir, outpath, "_", cov_gen1_name, "_ne.txt.done", sep="")
        if (file.exists(repl_ne_outpath_done)) {
            est_nes[repl] = scan(repl_ne_outpath)
        } else {
            # print("est_nes inputs:")
            # print(myTraj_repltemp[,traj_gen1_name])
            # print(myTraj_repltemp[,traj_gen2_name])
            # print(myCov_repltemp[,cov_gen1_name])
            # print(myCov_repltemp[,cov_gen2_name])
            # print(info$gen_levels[length(info$gen_levels)] - info$gen_levels[1])
            # print(pool1)
            # print(pool2)
            est_nes[repl] = estimateNe(
                p0=myTraj_repltemp[,traj_gen1_name], 
                pt=myTraj_repltemp[,traj_gen2_name], 
                cov0=myCov_repltemp[,cov_gen1_name], 
                covt=myCov_repltemp[,cov_gen2_name], 
                t=info$gen_levels[length(info$gen_levels)] - info$gen_levels[1],
                method=c("P.planII"),
                poolSize=c(pool1, pool2)
            )
            write(est_nes[repl], repl_ne_outpath, sep = "\t")
            write(est_nes[repl], repl_ne_outpath_done, sep = "\t")
        }
        # note: add options when ready: Ncensus=1000, poolSize=c(300, 300)
    }
    
    ne_outpath = paste(outdir, outpath, "_ne.txt", sep="")
    ne_outpath_done = paste(outdir, outpath, "_ne.txt.done", sep="")
    write(est_nes, ne_outpath, sep = "\t")
    write(est_nes, ne_outpath_done, sep = "\t")
    
    mean_ne = mean(est_nes)
    mean_ne_outpath = paste(outdir, outpath, "_mean_ne.txt", sep="")
    mean_ne_outpath_done = paste(outdir, outpath, "_mean_ne.txt.done", sep="")
    write(mean_ne, mean_ne_outpath, sep = "\t")
    write(mean_ne, mean_ne_outpath_done, sep = "\t")
    
    # out = est_full(mySync, mean_ne, info)
    out = est_full_save(mySync, mean_ne, info, outpath)
    # print(out)
    write.table(out, outpath, sep="\t", quote=FALSE, row.names=FALSE)
}

main()
