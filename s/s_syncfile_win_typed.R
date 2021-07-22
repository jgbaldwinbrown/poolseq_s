#!/usr/bin/env Rscript

# .onLoad <- function(libname, pkgname) {
#   backports::import(pkgname, c("deparse1"))
# }
# library("r-lib/backports")
library(utils)
library(backports)
deparse1 = getFromNamespace("deparse1", "backports")

suppressMessages(library(poolSeq))
suppressMessages(library(parallel))
library(typed)
# options(error = recover)
options(error = quote({dump.frames(to.file=TRUE); q()}))


global_chunksize = 1000

get_info <- function(infopath) {
    Data.frame() ? info_unstructured <- as.data.frame(fread(infopath, sep="\t", header=TRUE))
    # print(info_unstructured)
    List() ? info <- vector(mode = "list", length = 10)
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
    traj = af.traj(sync, chrom, pos, repl=info$repl_levels)
    est_p <- estimateSH(traj, Ne=round(Ne), t=info$gen_levels, h=0.5, simulate.p.value=TRUE)
    return(est_p)
}
estimateSH_one_win = function(sync, Ne, info, chrompos, one_start, winsize) {
    # check this
    # print("chrompos:")
    # print(chrompos)
    # print("one_start:")
    # print(one_start)
    span = one_start$index[1]:min(one_start$index[1] + winsize - 1, max(chrompos$index))
    chrom = chrompos$chrom[span]
    pos = chrompos$pos[span]
    # print("chrom:")
    # print(chrom)
    # print("pos:")
    # print(pos)
    traj = af.traj(sync, chrom, pos, repl=info$repl_levels)
    est_p <- estimateSH(traj, Ne=round(Ne), t=info$gen_levels, h=0.5, simulate.p.value=TRUE)
    return(est_p)
}

estimateSH_individual_loci <- function(sync, Ne, info) {
    chrom_pos = cbind(info$chrom, info$pos)
    chrom_pos_list = split(chrom_pos, seq(nrow(chrom_pos)))
    s_vals = unname(mclapply(chrom_pos_list, function(x) {estimateSH_one_locus(sync, Ne, info, x[1], x[2])}))
    return(list(chrom_pos, s_vals))
}

combine_est_sh_outputs <- function(full_output_list) {
    chrom_pos_list = mclapply(full_output_list, function(x) x[[1]])
    s_vals_list = unlist(mclapply(full_output_list, function(x) x[[2]]), recursive=FALSE)
    chrom_pos = do.call("rbind", chrom_pos_list)
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
            file.create(temppath_est_sh_done)
        } else {
            full_output_list[[j]] = readRDS(temppath_est_sh)
        }
        j = j + 1;
    }
    return(combine_est_sh_outputs(full_output_list))
}

estimateSH_wins <- function(sync, Ne, info, chrompos, start_chrompos_chunk, winsize) {
    # print("nrow(start_chrompos_chunk):")
    # print(nrow(start_chrompos_chunk))
    start_chrompos_chunk_list = split(start_chrompos_chunk, 1:nrow(start_chrompos_chunk))
    s_vals = unname(mclapply(start_chrompos_chunk_list, function(x) {estimateSH_one_win(sync, Ne, info, chrompos, x, winsize)}))
    # print("length(s_vals):")
    # print(length(s_vals))
    return(list(start_chrompos_chunk, s_vals))
}

# estimateSH_win_chunks <- function(sync, Ne, info, chrompos, all_start_chrompos, winsize) {
#     start_chrompos_list = split(all_start_chrompos, seq(nrow(all_start_chrompos)))
#     s_vals_chunks = vector(mode = "list", length = length(start_chrompos_list))
#     print("start_chrompos_list:")
#     print(start_chrompos_list)
#     for (i in 1:length(start_chrompos_list)) {
#         s_vals_chunks[i] = estimateSH_wins(sync, Ne, info, chrompos, start_chrompos_list[[i]], winsize)
#     }
#     s_vals = unname(unlist(s_vals_chunks))
#     return(list(all_start_chrompos, s_vals))
# }

estimateSH_win_savewrapper <- function(sync, Ne, info, outpath, winsize, winstep) {
    chrompos = data.frame(chrom = info$chrom, pos = as.numeric(info$pos), index = as.numeric(1:length(info$chrom)))
    chrompos = chrompos[order(chrompos[,"chrom"], chrompos["pos"]),]
    write.table(chrompos, "chrompos_temp.txt", sep="\t")
    all_starts = chrompos[seq(1,nrow(chrompos), winstep),]
    # print("all_starts:")
    # print(all_starts)
    #all_starts_chunked = split(all_starts, (as.numeric(rownames(all_starts))-1) %/% global_chunksize)
    all_starts_chunked = split(all_starts, (0:(nrow(all_starts)-1)) %/% global_chunksize)
    # print("all_starts_chunked:")
    # print(all_starts_chunked)
    full_output_list = vector(mode = "list", length = length(all_starts_chunked))
    j = 1
    for (i in 1:length(all_starts_chunked)) {
        temppath_est_sh = paste(outpath, "_tempdir/", outpath, "_est_sh_", as.character(i), ".RData", sep="")
        temppath_est_sh_done = paste(outpath, "_tempdir/", outpath, "_est_sh_", as.character(i), ".RData.done", sep="")
        if (! file.exists(temppath_est_sh_done)) {
            start_chunk = all_starts_chunked[[i]]
            # print("start_chunk:")
            # print(start_chunk)
            temp = estimateSH_wins(sync, Ne, info, chrompos, start_chunk, winsize)
            full_output_list[[j]] = temp
            saveRDS(temp, file = temppath_est_sh)
            file.create(temppath_est_sh_done)
        } else {
            full_output_list[[j]] = readRDS(temppath_est_sh)
        }
        j = j + 1;
    }
    out = combine_est_sh_outputs(full_output_list)
    # print("estimateSH_win_savewrapper length(out):")
    # print(length(out))
    # print("estimateSH_win_savewrapper nrow(out):")
    # print(nrow(out))
    return(out)
}

gets_from_ests <- function(ests) {
    # print("ests:")
    # print(ests)
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

est_full_save_win <- function(sync, Ne, info, outpath, winsize, winstep) {
    est_list <- estimateSH_win_savewrapper(sync, Ne, info, outpath, winsize, winstep)
    # print("est_list:")
    # print(est_list)
    # print("str(est_list):")
    # print(str(est_list))
    chrom_pos = est_list[[1]]
    est_all_p = est_list[[2]]
    s_vals = gets_from_ests(est_all_p)
    p_vals = getp_from_ests(est_all_p)
    # print(s_vals)
    # print(p_vals)
    out = as.data.frame(cbind(chrom_pos[,c("chrom", "pos")], s_vals, p_vals))
    # print(out)
    colnames(out) = c("chrom", "pos", "s", "p.value")
    # print("nrow(out)")
    # print(nrow(out))
    return(out)
}

est_full_save_repls_win <- function(sync, Ne, info, outpath, winsize, winstep) {
    out = vector(mode="list", length = length(info$repl_levels))
    for (i in 1:length(out)) {
        temp_info = info
        temp_info$repl_levels = info$repl_levels[i]
        outdir = paste(outpath, "_repl", as.character(temp_info$repl_levels), "_tempdir/", sep="")
        if (! dir.exists(outdir)) {
            dir.create(outdir)
        }

        repl_outpath = paste(outpath, "_repl", as.character(temp_info$repl_levels), sep="")
        est_list <- estimateSH_win_savewrapper(sync, Ne, temp_info, repl_outpath, winsize, winstep)
        # print("est_list:")
        # print(est_list)
        # print("str(est_list):")
        # print(str(est_list))
        chrom_pos = est_list[[1]]
        est_all_p = est_list[[2]]
        s_vals = gets_from_ests(est_all_p)
        p_vals = getp_from_ests(est_all_p)
        # print(s_vals)
        # print(p_vals)
        out[[i]] = as.data.frame(cbind(chrom_pos[,c("chrom", "pos")], s_vals, p_vals))
        # print(out)
        colnames(out[[i]]) = c("chrom", "pos", "s", "p.value")
        # print("nrow(out)")
        # print(nrow(out))
    }
    return(out)
}

est_full_save_repls <- function(sync, Ne, info, outpath) {
    out = vector(mode="list", length = length(info$repl_levels))
    for (i in 1:length(out)) {
        temp_info = info
        temp_info$repl_levels = info$repl_levels[i]
        # print(info)
        # print("full info:")
        # print(temp_info)
        outdir = paste(outpath, "_repl", as.character(temp_info$repl_levels), "_tempdir/", sep="")
        if (! dir.exists(outdir)) {
            dir.create(outdir)
        }
        est_list <- estimateSH_individual_loci_savewrapper(sync, Ne, temp_info, paste(outpath,"_repl",as.character(temp_info$repl_levels), sep=""))
        # print(est_list)
        chrom_pos = est_list[[1]]
        est_all_p = est_list[[2]]
        s_vals = gets_from_ests(est_all_p)
        p_vals = getp_from_ests(est_all_p)
        # print(s_vals)
        # print(p_vals)
        out[[i]] = as.data.frame(cbind(chrom_pos[,c("chrom", "pos")], s_vals, p_vals))
        # print(out)
        colnames(out[[i]]) = c("chrom", "pos", "s", "p.value")
    }
    return(out)
}

main = function() {
    
    args = commandArgs(trailingOnly = TRUE)
    syncpath = args[1]
    infopath = args[2]
    outpath = args[3]
    winsize = strtoi(args[4])
    winstep = strtoi(args[5])
    info = get_info(infopath)
    mySync <- read.sync(file=syncpath, gen=info$gen, repl=info$repl)
    info = update_info(info, mySync)
    
    outdir = paste(outpath, "_tempdir/", sep="")
    if (! dir.exists(outdir)) {
        dir.create(outdir)
    }
    
    est_nes = rep(NA, length(info$repl_levels))
    for (repl in info$repl_levels) {
        myTraj_repltemp = af.traj(mySync, info$chrom, info$pos, repl)
        myCov_repltemp = coverage(mySync, info$chrom, info$pos, repl=repl, gen=info$gen_levels)
        gen1 = info$gen_levels[1]
        gen2 = info$gen_levels[length(info$gen_levels)]
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
    
    out = est_full_save_win(mySync, mean_ne, info, outpath, winsize, winstep)
    out_repls = est_full_save_repls_win(mySync, mean_ne, info, outpath, winsize, winstep)
    write.table(out, outpath, sep="\t", quote=FALSE, row.names=FALSE)
    for (i in 1:length(out_repls)) {
        write.table(out_repls[[i]], paste(outpath, "_repl", info$repl_levels[i], sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    }
}

main()
