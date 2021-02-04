#!/bin/bash
set -e

export MC_CORES=1

exec ./s_syncfile.R \
    ex10.sync \
    ex10.info \
    ex10.tseries
