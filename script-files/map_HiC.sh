#!/usr/bin/bash

#SBATCH --cpus-per-task=38
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --time=8:00:00
#SBATCH --qos=short
#SBATCH --mem=65g
##SBATCH --constraint="c1"

set -u
###################################################################################################
#extract variables
VARI=$1
splitVARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$splitVARI"
echo $1 | tr ',' '\n'

TIME=$(date "+%s")
source ${SCRIPTdir}tools
LC_ALL=C

#determine available memory
MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/);{print X[5]}}')
MEM=$(($MEM - 5))
###################################################################################################
#setup-phase
TMPdirRAW=$TMPdir

locTMP=${TMPdir}prepareHiC/
mkdir -p $locTMP
echo $locTMP

###################################################################################################

cd $locTMP
bwa mem -A1 -B4  -E50 -L0 -t $SLURM_CPUS_PER_TASK ${locTMP}bwa-index ${locTMP}raw.${SLURM_ARRAY_TASK_ID}.fa | samtools view -Shb - > ${locTMP}mate_R${SLURM_ARRAY_TASK_ID}.bam