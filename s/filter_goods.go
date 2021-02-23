package main

import (
    "fmt"
    "os"
    "bufio"
    "strings"
    "strconv"
)

type Pos_key struct {
    chrom string
    pos int
}

type Entry []string

func get_pos_map_from_bed(path string) map[Pos_key]Entry {
    good_map := make(map[Pos_key]Entry)

    file, err := os.Open(path)
    if err != nil {
        panic(err)
    }

    scanner := bufio.NewScanner(file)
    for scanner.Scan() {
        split_line := strings.Split(strings.TrimSuffix(scanner.Text(), "\n"), "\t")
        chrom := split_line[0]
        pos, err := strconv.Atoi(split_line[2])
        if err != nil {
            panic(err)
        }
        good_map[Pos_key{chrom: chrom, pos: pos}] = split_line
    }
    return good_map
}

func print_goods(good_map map[Pos_key]Entry, inconn *os.File) {
    scanner := bufio.NewScanner(inconn)
    for scanner.Scan() {
        split_line := strings.Split(strings.TrimSuffix(scanner.Text(), "\n"), "\t")
        chrom := split_line[0]
        if chrom == "chrom" {
            fmt.Println(scanner.Text())
        } else {
            pos, err := strconv.Atoi(split_line[1])
            if err != nil {
                panic(err)
            }
            _, ok := good_map[Pos_key{chrom: chrom, pos: pos}]
            if ok {
                fmt.Println(scanner.Text())
            }
        }
    }
}

func main() {
    good_map := get_pos_map_from_bed(os.Args[1])
    print_goods(good_map, os.Stdin)
}
