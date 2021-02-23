package main

import (
    "fmt"
    "os"
    "bufio"
    "strings"
)

type Info_entry struct {
    path string
    breed string
    bitted string
    replicate string
}

func get_combine_info(inconn *os.File) []Info_entry {
    var info []Info_entry
    scanner := bufio.NewScanner(inconn)
    for scanner.Scan() {
        split_line := strings.Split(strings.TrimSuffix(scanner.Text(), "\n"), "\t")
        info = append(info, Info_entry{ path: split_line[0], breed: split_line[1], bitted: split_line[2], replicate: split_line[3]})
    }
    return info
}

func print_combo(info []Info_entry) {
    for file_number, entry := range info {
        file, err := os.Open(entry.path)
        if err != nil {
            panic(err)
        }
        scanner := bufio.NewScanner(file)
        line_number := 0
        for scanner.Scan() {
            line := strings.TrimSuffix(scanner.Text(), "\n")
            if line_number < 1 {
                if file_number < 1 {
                    fmt.Printf("%v\tbreed\tbitted\treplicate\tchrpos\tfull_treatment\n", line)
                }
            } else {
                split_line := strings.Split(line, "\t")
                fmt.Printf("%v\t%v\t%v\t%v\t%v:%v\t%v:%v\n", line, entry.breed, entry.bitted, entry.replicate, split_line[0], split_line[1], entry.breed, entry.bitted)
            }
            line_number++
        }
        cerr := file.Close()
        if cerr != nil {
            panic(cerr)
        }
    }
}

func main() {
    combine_info := get_combine_info(os.Stdin)
    print_combo(combine_info)
}
