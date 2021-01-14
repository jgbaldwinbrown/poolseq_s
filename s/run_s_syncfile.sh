#!/bin/bash
#SBATCH -t 72:00:00    #max:    72 hours (24 on ash)
#SBATCH -N 1          #format: count or min-max
#SBATCH -A owner-guest    #values: yandell, yandell-em (ember), ucgd-kp (kingspeak)
#SBATCH -p kingspeak-guest    #kingspeak, ucgd-kp, kingspeak-freecycle, kingspeak-guest
#SBATCH -J np_combo        #Job name

set -e

NUM_CORES="${SLURM_CPUS_ON_NODE}"
export MC_CORES="$NUM_CORES"

s_syncfile \
    SYNC_IN \
    INFO_IN \
> TIME_SERIES_OUT
touch TIME_SERIES_OUT.done
