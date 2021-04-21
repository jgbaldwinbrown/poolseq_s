#!/bin/bash
set -e

export MC_CORES=8

./s_syncfile_win.R \
    ex10.sync \
    ex10.info \
    ex10_win_3_2.tseries \
    3 \
    2

./s_syncfile_win.R \
    ex12.sync \
    ex12.info \
    ex12_win_3_2.tseries \
    3 \
    2

./s_syncfile_win.R \
    ex13.sync \
    ex13.info \
    ex13_win_3_2.tseries \
    3 \
    2

./s_syncfile_win.R \
    ex14.sync \
    ex14.info \
    ex14_win_3_2.tseries \
    3 \
    2

./combine_goods <ex10_win_goods.txt | \
tee combined_win_goods.txt | \
./tukey_s.py tukey_prefix_win
