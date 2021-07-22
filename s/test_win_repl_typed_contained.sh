#!/bin/bash
set -e

export MC_CORES=8

mkdir -p test_tukey_typed && cd test_tukey_typed && (
	cp ../ex10.sync .
	cp ../ex10.info .
	cp ../ex12.sync .
	cp ../ex12.info .
	cp ../ex13.sync .
	cp ../ex13.info .
	cp ../ex14.sync .
	cp ../ex14.info .
	cp ../s_syncfile_win_typed.R .
	cp ../ex10_win_goods.txt .
	cp ../tukey_s.py .
	cp ../combine_goods .
	
	./s_syncfile_win_typed.R \
		ex10.sync \
		ex10.info \
		ex10_win_3_2.tseries \
		3 \
		2
	
	./s_syncfile_win_typed.R \
		ex12.sync \
		ex12.info \
		ex12_win_3_2.tseries \
		3 \
		2
	
	./s_syncfile_win_typed.R \
		ex13.sync \
		ex13.info \
		ex13_win_3_2.tseries \
		3 \
		2
	
	./s_syncfile_win_typed.R \
		ex14.sync \
		ex14.info \
		ex14_win_3_2.tseries \
		3 \
		2
	
	./combine_goods <ex10_win_goods.txt | \
	tee combined_win_goods.txt | \
	./tukey_s.py tukey_prefix_win
)
