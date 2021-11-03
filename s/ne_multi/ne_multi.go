package main

import (
	"regexp"
	"runtime"
	"sync"
	"bufio"
	"github.com/montanaflynn/stats"
	"fmt"
	"os/exec"
	"os"
	"io/ioutil"
	"math/rand"
	"flag"
	"compress/gzip"
	"strconv"
	"strings"
)

var source string = `
#!/usr/bin/env Rscript

# .onLoad <- function(libname, pkgname) {
#   backports::import(pkgname, c("deparse1"))
# }
# library("r-lib/backports")
suppressMessages(library(utils))
suppressMessages(library(backports))
suppressMessages(library(data.table))
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
    List() ? info <- vector(mode = "list", length = 14)
    names(info) = c("chrom", "pos", "gen", "repl", "pool_size", "gen_levels", "repl_levels", "pool_size_levels", "chrom_levels", "pos_levels", "raw_chrom", "raw_pos", "raw_chrom_levels", "raw_pos_levels")
    info$gen = info_unstructured$gen
    info$repl = info_unstructured$repl
    info$pool_size = info_unstructured$pool_size
    info$gen_levels = sort(as.numeric(levels(factor(info_unstructured$gen))))
    info$repl_levels = sort(as.numeric(levels(factor(info_unstructured$repl))))
    info$pool_size_levels = sort(as.numeric(levels(factor(info_unstructured$pool_size))))
    # print(info)
    return(info)
}

update_info <- List() ? function(info= ? List(), sync, raw_sync= ? Data.frame()) {
    info$chrom = sync@alleles$chr
    info$pos = sync@alleles$pos
    info$chrom_levels = sort(levels(factor(info$chrom)))
    info$pos_levels = as.numeric(sort(levels(factor(info$pos))))

    info$raw_chrom = raw_sync[,1]
    info$raw_pos = raw_sync[,2]
    info$raw_chrom_levels = sort(levels(factor(info$raw_chrom)))
    info$raw_pos_levels = as.numeric(sort(levels(factor(info$raw_pos))))
    return(info)
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
    info = get_info(infopath)
    mySync <- read.sync(file=syncpath, gen=info$gen, repl=info$repl)
    myRawSync <- as.data.frame(fread(syncpath))
    info = update_info(info, mySync, myRawSync)
    
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
        repl_ne_outpath = paste(outpath, "_repl", as.character(repl), ".txt", sep="")
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
        # note: add options when ready: Ncensus=1000, poolSize=c(300, 300)
    }
    
    mean_ne = mean(est_nes)
    mean_ne_outpath = paste(outpath, ".txt", sep="")
    write(mean_ne, mean_ne_outpath, sep = "\t")
    
}

main()
`

var source_specify_gens string = `
#!/usr/bin/env Rscript

# .onLoad <- function(libname, pkgname) {
#   backports::import(pkgname, c("deparse1"))
# }
# library("r-lib/backports")
suppressMessages(library(utils))
suppressMessages(library(backports))
suppressMessages(library(data.table))
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
    List() ? info <- vector(mode = "list", length = 14)
    names(info) = c("chrom", "pos", "gen", "repl", "pool_size", "gen_levels", "repl_levels", "pool_size_levels", "chrom_levels", "pos_levels", "raw_chrom", "raw_pos", "raw_chrom_levels", "raw_pos_levels")
    info$gen = info_unstructured$gen
    info$repl = info_unstructured$repl
    info$pool_size = info_unstructured$pool_size
    info$gen_levels = sort(as.numeric(levels(factor(info_unstructured$gen))))
    info$repl_levels = sort(as.numeric(levels(factor(info_unstructured$repl))))
    info$pool_size_levels = sort(as.numeric(levels(factor(info_unstructured$pool_size))))
    # print(info)
    return(info)
}

update_info <- List() ? function(info= ? List(), sync, raw_sync= ? Data.frame()) {
    info$chrom = sync@alleles$chr
    info$pos = sync@alleles$pos
    info$chrom_levels = sort(levels(factor(info$chrom)))
    info$pos_levels = as.numeric(sort(levels(factor(info$pos))))

    info$raw_chrom = raw_sync[,1]
    info$raw_pos = raw_sync[,2]
    info$raw_chrom_levels = sort(levels(factor(info$raw_chrom)))
    info$raw_pos_levels = as.numeric(sort(levels(factor(info$raw_pos))))
    return(info)
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
    gen1 = as.numeric(args[3])
    gen2 = as.numeric(args[4])
    outpath = args[5]
    info = get_info(infopath)
    mySync <- read.sync(file=syncpath, gen=info$gen, repl=info$repl)
    myRawSync <- as.data.frame(fread(syncpath))
    info = update_info(info, mySync, myRawSync)
    
    outdir = paste(outpath, "_tempdir/", sep="")
    if (! dir.exists(outdir)) {
        dir.create(outdir)
    }
    
    est_nes = as.numeric(rep(NA, length(info$repl_levels)))
    for (repl in info$repl_levels) {
        myTraj_repltemp = af.traj(mySync, info$chrom, info$pos, repl)
        myCov_repltemp = coverage(mySync, info$chrom, info$pos, repl=repl, gen=info$gen_levels)
        # gen1 = info$gen_levels[1]
        # gen2 = info$gen_levels[length(info$gen_levels)]
        pool1 = info$pool_size[info$gen == gen1 & info$repl == repl]
        pool2 = info$pool_size[info$gen == gen2 & info$repl == repl]
        traj_gen1_name = paste("F", as.character(gen1), sep="")
        traj_gen2_name = paste("F", as.character(gen2), sep="")
        cov_gen1_name = paste("F", as.character(gen1), ".R", as.character(repl), ".cov", sep="")
        cov_gen2_name = paste("F", as.character(gen2), ".R", as.character(repl), ".cov", sep="")
        repl_ne_outpath = paste(outpath, "_repl", as.character(repl), ".txt", sep="")
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
        # note: add options when ready: Ncensus=1000, poolSize=c(300, 300)
    }
    
    mean_ne = mean(est_nes)
    mean_ne_outpath = paste(outpath, ".txt", sep="")
    write(mean_ne, mean_ne_outpath, sep = "\t")
    
}

main()
`

var rsource string = `
print(1)
`

type flag_set struct {
	randsource int
	data string
	perc_keep_string string
	perc_keep float64
	repls int
	info string
	oprefix string
	jobs int
}

type ne_t string

type ne_out_t struct {
	ne ne_t
	repl_nes []ne_t
	err error
}

type time_ne_t struct {
	time_0 string
	time_t string
	bio_repl string
	comp_repl int
	ne string
	err error
}

type time_ne_out_t struct {
	ne []time_ne_t
	repl_nes []time_ne_t
	err error
}

func get_flags() (flags flag_set) {
	flag.IntVar(&flags.randsource, "r", 0, "Random seed (integer).")
	flag.StringVar(&flags.data, "d", "", "Data path (gzipped sync file, must be a real file -- not a pipe).")
	flag.StringVar(&flags.perc_keep_string, "p", "", "Fraction of data to retain in each replicate.")
	flag.StringVar(&flags.info, "i", "", "Path of info file for this sync file.")
	flag.StringVar(&flags.oprefix, "o", "ne_multi_out", "Output prefix (default = ne_multi_out).")
	flag.IntVar(&flags.repls, "n", 1, "Number of replicates to perform.")
	flag.IntVar(&flags.jobs, "j", 1, "Maximum number of parallel jobs.")
	flag.Parse()
	if flags.data == "" || flags.perc_keep_string == "" || flags.info == "" {
		panic(fmt.Errorf("Error: missing argument.\n"))
	}

	var err error
	flags.perc_keep, err = strconv.ParseFloat(flags.perc_keep_string, 64)
	if err != nil { panic(err) }

	return flags
}

func subset_data(datapath string, perc_keep float64) (minipath string, err error) {
	minidataconn, err := ioutil.TempFile(".", "multi_ne.*.sync.gz")
	if err != nil { return minipath, err }
	defer minidataconn.Close()
	minipath = minidataconn.Name()

	minigzconn:= gzip.NewWriter(minidataconn)
	defer minigzconn.Close()

	dataconn, err := os.Open(datapath)
	if err != nil { return minipath, err }
	defer dataconn.Close()

	datagzconn, err := gzip.NewReader(dataconn)
	if err != nil { return minipath, err }
	defer datagzconn.Close()

	s := bufio.NewScanner(datagzconn)
	s.Buffer(make([]byte, 0), 1e12)

	for s.Scan() {
		rval := rand.Float64()
		if rval <= perc_keep {
			fmt.Fprintln(minigzconn, s.Text())
		}
	}

	return minipath, nil
}

func parse_ne(path string) (ne ne_t, err error) {
	err = nil
	conn, err := os.Open(path)
	if err != nil { return ne, err }
	defer conn.Close()
	s := bufio.NewScanner(conn)
	if s.Scan() {
		ne = ne_t(s.Text())
	} else {
		err = fmt.Errorf("Error: no Ne value.\n")
	}
	return ne, err
}

func parse_time_ne(path string) (ne time_ne_t, err error) {
	err = nil
	conn, err := os.Open(path)
	timere := regexp.MustCompile(`.*start([0-9]*).*end([0-9])`)
	if err != nil { return ne, err }
	defer conn.Close()
	s := bufio.NewScanner(conn)
	if s.Scan() {
		ne.ne = s.Text()
	} else {
		err = fmt.Errorf("Error: no Ne value.\n")
	}

	parsed := timere.FindStringSubmatch(path)
	ne.time_0 = parsed[1]
	ne.time_t = parsed[2]
	if err != nil { return ne, err }
	return ne, err
}

func parse_nes(paths ...string) (nes []ne_t, err error) {
	for _, path := range paths {
		ne, err := parse_ne(path)
		if err != nil { return nes, err }
		nes = append(nes, ne)
	}
	return nes, nil
}

func parse_time_nes(paths ...string) (time_nes []time_ne_t, err error) {
	for _, path := range paths {
		ne, err := parse_time_ne(path)
		if err != nil { return time_nes, err }
		time_nes = append(time_nes, ne)
	}
	return time_nes, nil
}

func calc_nes(datapath string, infopath string) (ne ne_t, repl_nes []ne_t, err error) {
	outdir, err := ioutil.TempDir(".", "ne_multi")
	if err != nil { return ne, repl_nes, err }
	defer os.RemoveAll(outdir)

	sourceconn, err := ioutil.TempFile(".", "ne_multi_source")
	if err != nil { return ne, repl_nes, err }
	sourcepath := sourceconn.Name()
	defer os.Remove(sourcepath)
	fmt.Fprintln(sourceconn, source)
	sourceconn.Close()

	command := exec.Command("Rscript", sourcepath, datapath, infopath, outdir + "/ne")
	err = command.Run()
	if err != nil { return ne, repl_nes, err }

	ne, err = parse_ne(outdir + "/ne.txt")
	if err != nil { return ne, repl_nes, err }
	repl_nes, err = parse_nes(
		outdir + "/ne_repl1.txt",
		outdir + "/ne_repl2.txt",
		outdir + "/ne_repl3.txt",
		outdir + "/ne_repl4.txt",
	)
	if err != nil { return ne, repl_nes, err }

	return ne, repl_nes, nil
}

func calc_time_nes(datapath string, infopath string) (time_ne []time_ne_t, repl_time_nes []time_ne_t, err error) {
	outdir, err := ioutil.TempDir(".", "ne_multi")
	if err != nil { return time_ne, repl_time_nes, err }
	defer os.RemoveAll(outdir)

	sourceconn, err := ioutil.TempFile(".", "ne_multi_time_source")
	if err != nil { return time_ne, repl_time_nes, err }
	sourcepath := sourceconn.Name()
	defer os.Remove(sourcepath)
	fmt.Fprintln(sourceconn, source)
	sourceconn.Close()

	times := []int{0, 6, 12, 18, 24, 30, 36, 42, 48, 54}
	var replpaths []string
	var meanpaths []string

	for time_t_i :=1; time_t_i < len(times); time_t_i++ {
		time_0 := times[time_t_i-1]
		time_t := times[time_t_i]
		outpref := fmt.Sprintf("%s/ne_start%v_end%v", outdir, time_0, time_t)
		command := exec.Command("Rscript", sourcepath, datapath, infopath, fmt.Sprintf("%v", time_0), fmt.Sprintf("%v", time_t), outdir + "/ne")
		err = command.Run()
		if err != nil { return time_ne, repl_time_nes, err }
		meanpaths = append(meanpaths, outpref + ".txt")
		for i:=1; i<5; i++ {
			replpath := fmt.Sprintf("%s_repl%v.txt", outpref, i)
			replpaths = append(replpaths, replpath)
		}
	}

	time_ne, err = parse_time_nes(meanpaths...)
	if err != nil { return time_ne, repl_time_nes, err }
	repl_time_nes, err = parse_time_nes(replpaths...)
	if err != nil { return time_ne, repl_time_nes, err }

	return time_ne, repl_time_nes, nil
}

func get_ne(datapath string, perc_keep float64, infopath string) (ne ne_t, repl_nes []ne_t, err error) {
	minipath, err := subset_data(datapath, perc_keep)
	defer os.Remove(minipath)
	if err != nil { return ne, repl_nes, err }

	ne, repl_nes, err = calc_nes(minipath, infopath)

	return ne, repl_nes, err
}

func get_time_ne(datapath string, perc_keep float64, infopath string) (time_ne []time_ne_t, repl_time_nes []time_ne_t, err error) {
	minipath, err := subset_data(datapath, perc_keep)
	defer os.Remove(minipath)
	if err != nil { return time_ne, repl_time_nes, err }

	time_ne, repl_time_nes, err = calc_time_nes(minipath, infopath)

	return time_ne, repl_time_nes, err
}

func get_ne_parallel(datapath string, perc_keep float64, infopath string, ne_out *ne_out_t, parent_wg *sync.WaitGroup) {
	defer parent_wg.Done()
	ne_out.ne, ne_out.repl_nes, ne_out.err = get_ne(datapath, perc_keep, infopath)
}

func get_time_ne_parallel(datapath string, perc_keep float64, infopath string, time_ne_out *time_ne_out_t, parent_wg *sync.WaitGroup) {
	defer parent_wg.Done()
	time_ne_out.ne, time_ne_out.repl_nes, time_ne_out.err = get_time_ne(datapath, perc_keep, infopath)
}

func split_outs(ne_outs []ne_out_t) (nes []ne_t, repl_nes [][]ne_t, err error) {
	for _, out := range ne_outs {
		nes = append(nes, out.ne)
		repl_nes = append(repl_nes, out.repl_nes)
		if out.err != nil { return nes, repl_nes, err }
	}
	return nes, repl_nes, nil
}

func split_time_ne_outs(ne_outs []time_ne_out_t) (nes []time_ne_t, repl_nes []time_ne_t, err error) {
	for _, out := range ne_outs {
		nes = append(nes, out.ne...)
		repl_nes = append(repl_nes, out.repl_nes...)
		if out.err != nil { return nes, repl_nes, err }
	}
	return nes, repl_nes, nil
}

func get_multi_nes(datapath string, randsource int, perc_keep float64, repls int, infopath string) (nes []ne_t, repl_nes [][]ne_t, err error) {
	rand.New(rand.NewSource(int64(randsource)))
	var wg sync.WaitGroup
	ne_outs := make([]ne_out_t, repls)
	for i:=0; i<repls; i++ {
		wg.Add(1)
		go get_ne_parallel(datapath, perc_keep, infopath, &ne_outs[i], &wg)
	}
	wg.Wait()
	nes, repl_nes, err = split_outs(ne_outs)
	return nes, repl_nes, err
}

func get_multi_time_nes(datapath string, randsource int, perc_keep float64, comp_repls int, infopath string) (time_nes []time_ne_t, repl_time_nes []time_ne_t, err error) {
	rand.New(rand.NewSource(int64(randsource)))
	var wg sync.WaitGroup
	time_ne_outs := make([]time_ne_out_t, comp_repls)
	for i:=0; i<comp_repls; i++ {
		wg.Add(1)
		go get_time_ne_parallel(datapath, perc_keep, infopath, &time_ne_outs[i], &wg)
	}
	wg.Wait()
	time_nes, repl_time_nes, err = split_time_ne_outs(time_ne_outs)
	return time_nes, repl_time_nes, err
}

func print_nes(nes []ne_t, path string) (err error) {
	conn, err := os.Create(path)
	if err != nil { return err }
	defer conn.Close()
	for _, ne := range nes {
		fmt.Fprintln(conn, ne)
	}
	return nil
}

func print_ne(ne ne_t, path string) (err error) {
	return print_nes([]ne_t{ne}, path)
}

func print_time_nes(nes []time_ne_t, path string) (err error) {
	conn, err := os.Create(path)
	if err != nil { return err }
	defer conn.Close()
	for _, ne := range nes {
		fmt.Fprintf(conn, "%v\t%v\t%v\t%v\t%v\n", ne.ne, ne.time_0, ne.time_t, ne.bio_repl, ne.comp_repl)
	}
	return nil
}

func floatify(nes []ne_t) (floats []float64) {
	for _, ne := range nes {
		float, err := strconv.ParseFloat(string(ne), 64)
		if err == nil {
			floats = append(floats, float)
		}
	}
	return floats
}

func print_repl_nes(nes [][]ne_t, path string) (err error) {
	conn, err := os.Create(path)
	if err != nil { return err }
	defer conn.Close()
	for _, ne_set := range nes {
		ne_set_strings := []string{}
		for _, ne := range ne_set {ne_set_strings = append(ne_set_strings, string(ne))}
		outstring := strings.Join(ne_set_strings, "\t")
		fmt.Fprintln(conn, outstring)
	}
	return nil
}

func floatify_repls( nes [][]ne_t) (all_floats [][]float64) {
	for _, ne_set := range nes {
		all_floats = append(all_floats, floatify(ne_set))
	}
	return all_floats
}

func mean_repls(nes [][]float64) (means []float64) {
	if len(nes) < 1 {
		return means
	}
	for i:=0; i<len(nes[0]); i++ {
		ne_set := []float64{}
		for _, nerow := range nes {
			if len(nerow) > i {
				ne_set = append(ne_set, nerow[i])
			}
		}
		mean, err := stats.Mean(ne_set)
		if err == nil {
			means = append(means, mean)
		}
	}
	return means
}

func mean_time_nes(time_nes []time_ne_t) (means []time_ne_t, err error) {
	times := []string{"0", "6", "12", "18", "24", "30", "36", "42", "48", "54"}
	if len(time_nes) < 1 {
		return means, fmt.Errorf("Error: length less than 1.\n")
	}
	for time_i, time := range times[:len(times)-1] {
		ne_set := []float64{}
		for _, ne := range time_nes {
			if ne.time_0 == time {
				ne_float, err := strconv.ParseFloat(ne.ne, 64)
				if err != nil { return means, err }
				ne_set = append(ne_set, ne_float)
			}
		}
		mean, err := stats.Mean(ne_set)
		if err != nil { return means, err }
		means = append(means, time_ne_t{ne: fmt.Sprintf("%v", mean), time_0: time, time_t : times[time_i], bio_repl: "-1", comp_repl: -1})
	}
	return means, nil
}

func mean_repl_time_nes(time_nes []time_ne_t) (means []time_ne_t, err error) {
	times := []string{"0", "6", "12", "18", "24", "30", "36", "42", "48", "54"}
	if len(time_nes) < 1 {
		return means, fmt.Errorf("Error: time_nes too short.")
	}
	for time_i, time := range times[:len(times)-1] {
		for _, bio_repl := range []string{"1", "2", "3", "4"} {
			ne_set := []float64{}
			for _, ne := range time_nes {
				if ne.time_0 == time && ne.bio_repl == bio_repl {
					ne_float, err := strconv.ParseFloat(ne.ne, 64)
					if err != nil {return means, err}
					ne_set = append(ne_set, ne_float)
				}
			}
			mean, err := stats.Mean(ne_set)
			if err != nil { return means, err }
			means = append(means, time_ne_t{ne: fmt.Sprintf("%v", mean), time_0: time, time_t : times[time_i], bio_repl: "-1", comp_repl: -1})
		}
	}
	return means, nil
}

func print_floats(floats []float64, path string) (err error) {
	conn, err := os.Create(path)
	if err != nil { return err }
	defer conn.Close()
	for _, float := range floats {
		fmt.Fprintln(conn, float)
	}
	return nil
}

func main() {
	flags := get_flags()

	native_jobs := runtime.GOMAXPROCS(-1)
	if flags.jobs < native_jobs {
		runtime.GOMAXPROCS(flags.jobs)
	}

	nes, repl_nes, err := get_multi_nes(flags.data, flags.randsource, flags.perc_keep, flags.repls, flags.info)
	if err != nil { panic(err) }

	nes_path := flags.oprefix + "_nes.txt"
	err = print_nes(nes, nes_path)
	if err != nil { panic(err) }

	ne_path := flags.oprefix + "_mean_ne.txt"
	mean, err := stats.Mean(floatify(nes))
	if err != nil { panic(err) }
	err = print_floats([]float64{mean}, ne_path)
	if err != nil { panic(err) }

	repl_nes_path := flags.oprefix + "_repl_nes.txt"
	err = print_repl_nes(repl_nes, repl_nes_path)
	if err != nil { panic(err) }
	repl_mean_ne_path := flags.oprefix + "_mean_repl_ne.txt"
	err = print_floats(mean_repls(floatify_repls(repl_nes)), repl_mean_ne_path)
	if err != nil { panic(err) }

	time_nes, repl_time_nes, err := get_multi_time_nes(flags.data, flags.randsource, flags.perc_keep, flags.repls, flags.info)
	if err != nil { panic(err) }

	time_nes_path := flags.oprefix + "_time_nes.txt"
	err = print_time_nes(time_nes, time_nes_path)
	if err != nil { panic(err) }

	time_ne_path := flags.oprefix + "_mean_time_ne.txt"
	time_means, err := mean_time_nes(time_nes)
	if err != nil { panic(err) }
	err = print_time_nes(time_means, time_ne_path)
	if err != nil { panic(err) }

	repl_time_nes_path := flags.oprefix + "_repl_time_nes.txt"
	err = print_time_nes(repl_time_nes, repl_time_nes_path)
	if err != nil { panic(err) }
	repl_mean_time_ne_path := flags.oprefix + "_mean_repl_time_ne.txt"
	mean_repl_time_nes_val, err := mean_repl_time_nes(repl_time_nes)
	if err != nil { panic(err) }
	err = print_time_nes(mean_repl_time_nes_val, repl_mean_time_ne_path)
	if err != nil { panic(err) }
}


/*
    for i := 0; i < len(move_vecs); i++ {
        new_state := s.update(move_vecs[i])
        wg.Add(1)
        go new_state.score_parallel(&wg, depth - 1, planned_moves, &return_scores[i])
    }

    wg.Wait()
*/
