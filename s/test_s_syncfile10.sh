#!/bin/bash
set -e

export MC_CORES=1

./s_syncfile.R \
    ex10.sync \
    ex10.info \
    ex10.tseries

./s_syncfile.R \
    ex12.sync \
    ex12.info \
    ex12.tseries

./s_syncfile.R \
    ex13.sync \
    ex13.info \
    ex13.tseries

./s_syncfile.R \
    ex14.sync \
    ex14.info \
    ex14.tseries
