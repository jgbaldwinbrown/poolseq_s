package main

import (
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
}

type ne_t string

func printexample() {
	sourceconn, err := ioutil.TempFile(".", "simple")
	sourcepath := sourceconn.Name()
	if err != nil { panic(err) }
	defer os.Remove(sourcepath)

	fmt.Fprintln(sourceconn, rsource)
	sourceconn.Close()

	cmd_tmp := exec.Command("cat", sourcepath)
	cmd_tmp.Stdout = os.Stdout
	cmd_tmp.Stderr = os.Stderr
	err = cmd_tmp.Run()

	cmd := exec.Command("Rscript", sourcepath)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	err = cmd.Run()
	if err != nil { panic(err) }
}

func get_flags() (flags flag_set) {
	flag.IntVar(&flags.randsource, "r", 0, "Random seed (integer).")
	flag.StringVar(&flags.data, "d", "", "Data path (gzipped sync file, must be a real file -- not a pipe).")
	flag.StringVar(&flags.perc_keep_string, "p", "", "Fraction of data to retain in each replicate.")
	flag.StringVar(&flags.info, "i", "", "Path of info file for this sync file.")
	flag.StringVar(&flags.oprefix, "o", "ne_multi_out", "Output prefix (default = ne_multi_out).")
	flag.IntVar(&flags.repls, "n", 1, "Number of replicates to perform.")
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

func parse_nes(paths ...string) (nes []ne_t, err error) {
	for _, path := range paths {
		ne, err := parse_ne(path)
		if err != nil { return nes, err }
		nes = append(nes, ne)
	}
	return nes, nil
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

func get_ne(datapath string, perc_keep float64, infopath string) (ne ne_t, repl_nes []ne_t, err error) {
	minipath, err := subset_data(datapath, perc_keep)
	defer os.Remove(minipath)
	if err != nil { return ne, repl_nes, err }

	ne, repl_nes, err = calc_nes(minipath, infopath)

	return ne, repl_nes, err
}

func get_multi_nes(datapath string, randsource int, perc_keep float64, repls int, infopath string) (nes []ne_t, repl_nes [][]ne_t, err error) {
	rand.New(rand.NewSource(int64(randsource)))
	for i:=0; i<repls; i++ {
		ne, repl_ne, err := get_ne(datapath, perc_keep, infopath)
		if err != nil { return nes, repl_nes, err }
		nes = append(nes, ne)
		repl_nes = append(repl_nes, repl_ne)
	}
	return nes, repl_nes, nil
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
}
