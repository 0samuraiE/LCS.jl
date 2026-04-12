#!/bin/bash
#$ -cwd
#$ -N N1500
#$ -o logs/$JOB_NAME/$JOB_ID.out
#$ -e logs/$JOB_NAME/$JOB_ID.err

#$ -l h_rt=24:00:00
#$ -l gpu_1=8

set -euo pipefail

mkdir -p logs/$JOB_NAME

module purge
module load openmpi/5.0.2-gcc
module load hdf5-parallel/1.14.3/gcc11.4.1
module load cuda/12.8.0
module load cudnn/9.8.0

export PATH=/gs/fs/tga-onilab/tominaga/.juliaup/bin${PATH:+:${PATH}}
export JULIA_DEPOT_PATH=/gs/fs/tga-onilab/tominaga/.julia
export UCX_WARN_UNUSED_ENV_VARS=n
export UCX_ERROR_SIGNALS=SIGILL,SIGBUS,SIGFPE

LCSCUDA="LCSCUDA"

CFG="cfgs/$(echo $JOB_NAME | tr '%' '/').lcs-yaml"

julia --project="$LCSCUDA" -e 'using Pkg; Pkg.precompile()'

mpirun \
    -n 8 \
    -npernode 1 \
    -x LD_LIBRARY_PATH \
    -x UCX_WARN_UNUSED_ENV_VARS \
    -x UCX_ERROR_SIGNALS \
    -x LCS_LOG_LEVEL=PROFILE \
    -x LCS_RESUME=1 \
    julia --project="$LCSCUDA" "$LCSCUDA/simulate.jl" "${CFG}"
