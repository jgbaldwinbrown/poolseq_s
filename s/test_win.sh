#!/bin/bash
set -e

export MC_CORES=8

./s_syncfile_win.R \
    ex15.sync \
    ex15.info \
    ex15_win_3_2.tseries \
    3 \
    2
