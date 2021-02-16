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

get_partial_est_full_save <- function(sync, Ne, info, outpath) {
    est_list <- estimateSH_individual_loci_savewrapper(sync, Ne, info, outpath)
    chrom_pos = est_list[[1]]
    est_all_p = est_list[[2]]
    s_vals = gets_from_ests(est_all_p)
    p_vals = getp_from_ests(est_all_p)
    out = as.data.frame(cbind(chrom_pos, s_vals, p_vals))
    colnames(out) = c("chrom", "pos", "s", "p.value")
    return(out)
}

get_partial_est_full_save_repls <- function(sync, Ne, info, outpath) {
    out = vector(mode="list", length = length(info$repl_levels))
    for (i in 1:length(out)) {
        temp_info = info
        temp_info$repl_levels = info$repl_levels[i]
        outdir = paste(outpath, "_repl", as.character(temp_info$repl_levels), "_tempdir/", sep="")
        if (! dir.exists(outdir)) {
            dir.create(outdir)
        }
        est_list <- estimateSH_individual_loci_savewrapper(sync, Ne, temp_info, paste(outpath,"_repl",as.character(temp_info$repl_levels), sep=""))
        chrom_pos = est_list[[1]]
        est_all_p = est_list[[2]]
        s_vals = gets_from_ests(est_all_p)
        p_vals = getp_from_ests(est_all_p)
        out[[i]] = as.data.frame(cbind(chrom_pos, s_vals, p_vals))
        colnames(out[[i]]) = c("chrom", "pos", "s", "p.value")
    }
    return(out)
}

main = function() {
    
    args = commandArgs(trailingOnly = TRUE)
    syncpath = args[1]
    infopath = args[2]
    outpath = args[3]
    table_outpath = paste(outpath, "_partial", sep="")
    info = get_info(infopath)
    mySync <- read.sync(file=syncpath, gen=info$gen, repl=info$repl)
    info = update_info(info, mySync)
    
    outdir = paste(outpath, "_tempdir/", sep="")
    if (! dir.exists(outdir)) {
        dir.create(outdir)
    }
    
    est_nes = rep(NA, length(info$repl_levels))
    # out = est_full(mySync, mean_ne, info)
    out = get_partial_est_full_save(mySync, mean_ne, info, outpath)
    # out_repls = get_partial_est_full_save_repls(mySync, mean_ne, info, outpath)
    write.table(out, table_outpath, sep="\t", quote=FALSE, row.names=FALSE)
    # for (i in 1:length(out_repls)) {
    #     write.table(out_repls[[i]], paste(table_outpath, "_repl", info$repl_levels[i], sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    # }
}

main()
