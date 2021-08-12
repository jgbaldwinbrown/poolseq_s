#!/usr/bin/env Rscript

# .onLoad <- function(libname, pkgname) {
#   backports::import(pkgname, c("deparse1"))
# }
# library("r-lib/backports")
suppressMessages(library(utils))
suppressMessages(library(backports))
deparse1 = getFromNamespace("deparse1", "backports")

suppressMessages(library(poolSeq))
suppressMessages(library(parallel))
suppressMessages(library(typed))
# options(error = recover)
options(error = quote({dump.frames(to.file=TRUE); q()}))


global_chunksize = 1000

get_info <- List() ? function(infopath= ? Character()) {
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

update_info <- List() ? function(info= ? List(), sync) {
    info$chrom = sync@alleles$chr
    info$pos = sync@alleles$pos
    info$chrom_levels = sort(levels(factor(info$chrom)))
    info$pos_levels = as.numeric(sort(levels(factor(info$pos))))
    return(info)
}

estimateSH_one_win = function(sync, Ne= ? Double(), info= ? List(), chrompos= ? Data.frame(), one_start= ? Data.frame(), winsize= ? Integer()) {
    # check this
    # print("chrompos:")
    # print(chrompos)
    # print("one_start:")
    # print(one_start)
    Integer() ? span
    Double() ? pos
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

combine_est_sh_outputs <- List(2) ? function(full_output_list= ? List()) {
    List() ? chrom_pos_list
    List() ? s_vals_list
    Data.frame() ? chrom_pos
    chrom_pos_list = mclapply(full_output_list, function(x) x[[1]])
    s_vals_list = unlist(mclapply(full_output_list, function(x) x[[2]]), recursive=FALSE)
    chrom_pos = do.call("rbind", chrom_pos_list)
    return(list(chrom_pos, s_vals_list))
}

estimateSH_wins <- List(2) ? function(sync, Ne= ? Double(), info= ? List(), chrompos= ? Data.frame(), start_chrompos_chunk= ? Data.frame(), winsize= ? Integer()) {
    # print("nrow(start_chrompos_chunk):")
    # print(nrow(start_chrompos_chunk))
    List() ? start_chrompos_chunk_list
    start_chrompos_chunk_list = split(start_chrompos_chunk, 1:nrow(start_chrompos_chunk))
    List() ? s_vals
    s_vals = unname(mclapply(start_chrompos_chunk_list, function(x) {estimateSH_one_win(sync, Ne, info, chrompos, x, winsize)}))
    # print("length(s_vals):")
    # print(length(s_vals))
    return(list(start_chrompos_chunk, s_vals))
}

estimateSH_win_savewrapper <- List(2) ? function(sync, Ne= ? Double(), info= ? List(), outpath= ? Character(1), winsize= ? Integer(1), winstep= ? Integer(1)) {
    Data.frame() ? chrompos
    Data.frame() ? all_starts
    List() ? all_starts
    List() ? full_output_list
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
	Character(1) ? temppath_est_sh
	Character(1) ? temppath_est_sh_done
        temppath_est_sh = paste(outpath, "_tempdir/", outpath, "_est_sh_", as.character(i), ".RData", sep="")
        temppath_est_sh_done = paste(outpath, "_tempdir/", outpath, "_est_sh_", as.character(i), ".RData.done", sep="")
        if (! file.exists(temppath_est_sh_done)) {
            Data.frame() ? start_chunk
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
    List(2) ? out
    out = combine_est_sh_outputs(full_output_list)
    # print("estimateSH_win_savewrapper length(out):")
    # print(length(out))
    # print("estimateSH_win_savewrapper nrow(out):")
    # print(nrow(out))
    return(out)
}

gets_from_ests <- Double() ? function(ests= ? List()) {
    # print("ests:")
    # print(ests)
    List() ? s_vals
    # print("ests:")
    # print(ests)
    s_vals = mclapply(ests, function(x) {x$s})
    # print("s_vals:")
    # print(s_vals)
    return(unlist(s_vals))
}

getp_from_ests <- Double() ? function(ests= ? List()) {
    List() ? p_vals
    p_vals = mclapply(ests, function(x) {x$p.value})
    return(unlist(p_vals))
}

est_full_save_win <- Data.frame() ? function(sync, Ne= ? Double(), info= ? List(), outpath= ? Character(1), winsize= ? Integer(1), winstep= ? Integer(1)) {
    List(2) ? est_list
    est_list <- estimateSH_win_savewrapper(sync, Ne, info, outpath, winsize, winstep)
    List() ? chrom_pos
    List() ? est_all_p
    Double() ? s_vals
    Double() ? p_vals
    Data.frame() ? out
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

est_full_save_repls_win <- List() ? function(sync, Ne= ? Double(), info= ? List(), outpath= ? Character(1), winsize= ? Integer(1), winstep= ? Integer(1)) {
    List() ? out
    List() ? temp_info
    Character(1) ? outdir
    Character(1) ? repl_outpath
    List(2) ? est_list
    List() ? chrom_pos
    List() ? est_all_p
    Double() ? s_vals
    Double() ? p_vals
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

main = function() {
    
    Character() ? args
    Character(1) ? syncpath
    Character(1) ? infopath
    Character(1) ? outpath
    Integer(1) ? winsize
    Integer(1) ? sinstep
    List() ? info
    Character(1) ? outdir
    Double() ? est_nes
    Double() ? myTraj_repltemp
    Double () ? myCov_repltemp
    Double(1) ? gen1
    Double(1) ? gen2
    Integer(1) ? pool1
    Integer(1) ? pool2
    Character(1) ? traj_gen1_name
    Character(1) ? traj_gen2_name
    Character(1) ? cov_gen1_name
    Character(1) ? cov_gen2_name
    Character(1) ? repl_ne_outpath
    Character(1) ? repl_ne_outpath_done
    Character(1) ? ne_outpath
    Character(1) ? ne_outpath_done
    Double(1) ? mean_ne
    Character(1) ? mean_ne_outpath
    Character(1) ? mean_ne_outpath_done
    Data.frame() ? out
    List() ? out_repls

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
    
    est_nes = as.numeric(rep(NA, length(info$repl_levels)))
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
