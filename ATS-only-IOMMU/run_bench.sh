#!/bin/bash

ARGC=$#

# Benchmark name
BENCH_PREFIX="TDLCA"
BENCH_SUFFIX="td"

if [ $ARGC != 1 ]; then
    echo "Invalid argument number! ./run_bench.sh [benchmark name]"
    exit
fi

BENCH_NAME=$1
BENCH_LIST=(
#### Medical Imaging #####
    'Denoise'
    'Segmentation'
#### Commercial ##########
    'BlackScholes'
    'StreamCluster'
    'Swaptions'
#### Computer vision #####
    'Disparity_Map'
#### Navigation ##########
    'Robot_Localization'
    'EKF_SLAM'
)

FOUND="0"
for((i=0; i < ${#BENCH_LIST[@]}; i++)) do
    if [ ${BENCH_LIST[$i]} == $BENCH_NAME ]; then
        FOUND="1"
	break
    fi
done

if [ $FOUND == "0" ]; then
    echo "Invalid benchmarks! [Denoise, Segmentation, BlackScholes, StreamCluster, Swaptions, Disparity_Map, Robot_Localization, EKF_SLAM]"
    exit
fi

# Set up paths
RUN_BIN=$M5_PATH'/build/X86/gem5.opt'
DIR_OUT=$M5_PATH'/result'
DIR_CFG=$M5_PATH'/configs'
DIR_FSBOOT=$DIR_CFG'/boot'
DIR_CKPT=$M5_PATH'/ckpt-1core'

# Configuration setup
DEBUG_FLAG="--debug-flags=PageTableWalker"
CFG_OPTIONS="${DIR_CFG}/example/fs_tlb.py --checkpoint-dir=${DIR_CKPT} \
    --restore-with-cpu=detailed -r 1 -n 1 --l2_size=64kB --num-l2caches=32 \
    --mem-size=2GB --num-dirs=4 --ruby --lcacc --garnet=fixed --topology=Mesh \
    --mesh-rows=4 --num_accinstances=8 --host_ptw_latency=1 --lcacc_tlb_size=32 \
    --lcacc_tlb_mshr=0"

# Prepare and run command
FILE_BOOTSCRIPT=$DIR_FSBOOT'/'$BENCH_NAME'.'$BENCH_SUFFIX'.rcS'
DIR_OUT_BENCH=$DIR_OUT'/'$BENCH_PREFIX'_'$BENCH_NAME
FILE_OUT=$DIR_OUT_BENCH'/result.txt'

mkdir -p $DIR_OUT_BENCH

echo "(time -p ${RUN_BIN} --outdir=${DIR_OUT_BENCH} \
    ${CFG_OPTIONS} --acc_type=${BENCH_NAME} --work-end-exit-count=1 \
    --script=${FILE_BOOTSCRIPT}) >& ${FILE_OUT} &"
(time -p ${RUN_BIN} --outdir=${DIR_OUT_BENCH} \
    ${CFG_OPTIONS} --acc_type=${BENCH_NAME} --work-end-exit-count=1 \
    --script=${FILE_BOOTSCRIPT}) >& ${FILE_OUT} &
