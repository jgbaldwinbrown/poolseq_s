#!/bin/bash
set -e

export MC_CORES=8

exec ./s_syncfile.R \
    ex8.sync \
    ex8.info
