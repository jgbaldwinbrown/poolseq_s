package main

import (
	"github.com/jgbaldwinbrown/fasttsv"
	"io"
	"os"
	"strconv"
)

func int_entries(strings []string, ints []int) ([]int, error) {
	ints = ints[:0]
	for _,s := range strings {
		var err error
		i, err := strconv.Atoi(s)
		if err != nil {return ints, err}
		ints = append(ints, i)
	}
	return ints, nil
}

func bad_values(ints []int) bool {
	total_nonzero := 0
	for _,i := range ints {
		if i>0 {
			total_nonzero++
		}
	}
	return total_nonzero < 2
}

func good_line(l []string) (ok bool, bads []int) {
	var entry_str_arr [5]string
	var entry_int_arr [5]int
	entry_str := entry_str_arr[:0]
	entry_int := entry_int_arr[:0]
	total_nonzero := 0
	for i, e := range l[3:] {
		entry_str = fasttsv.ManualSplit(e, ':', entry_str)
		entry_int, err := int_entries(entry_str, entry_int)
		if err == nil && !bad_values(entry_int[:4]) {
			total_nonzero++
		} else {
			bads = append(bads, i)
		}
	}
	ok =  total_nonzero >= 2
	return ok, bads
}

func print_bad_lines(r io.Reader, W io.Writer) {
	s := fasttsv.NewScanner(r)
	w := fasttsv.NewWriter(W)
	for s.Scan() {
		ok, _ := good_line(s.Line())
		if ! ok {
			w.Write(s.Line())
			// fmt.Println(s.Line(), bads)
		}
	}
	w.Flush()
}

func main() {
	print_bad_lines(os.Stdin, os.Stdout)
}
