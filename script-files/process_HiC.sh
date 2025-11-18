#!/usr/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --time=24:00:00
#SBATCH --qos=medium
#SBATCH --mem=140g
##SBATCH --partition=m
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

#set THREADS
NUMEXPR_MAX_THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

#filter bam-files
samtools view -h ${locTMP}mate_R1.bam | filter_bam | samtools view -@ $SLURM_CPUS_PER_TASK -Sb - >${locTMP}mate_R1.filtered.bam &
samtools view -h ${locTMP}mate_R2.bam | filter_bam | samtools view -@ $SLURM_CPUS_PER_TASK -Sb - >${locTMP}mate_R2.filtered.bam &
wait

#combine bam-files
samtools faidx ${assemblyFASTA}
combine_bam ${locTMP}mate_R1.filtered.bam ${locTMP}mate_R2.filtered.bam samtools 10 | samtools view -bS -t ${assemblyFASTA}.fai - | samtools sort -@ $THREADS -o ${locTMP}combined.bam  -


#add read-group
AddOrReplaceReadGroups -INPUT ${locTMP}combined.bam -OUTPUT ${locTMP}combined.paired.bam -ID combined.paired -LB combined.paired -SM Scaffolding -PL ILLUMINA -PU none

#mark duplicates
MarkDuplicates -INPUT ${locTMP}combined.paired.bam -OUTPUT ${locTMP}combined.polished.bam -METRICS_FILE ${locTMP}metrics.txt -TMP_DIR $locTMP -ASSUME_SORTED TRUE -VALIDATION_STRINGENCY LENIENT -REMOVE_DUPLICATES TRUE

#index and create statistics
samtools index ${locTMP}combined.polished.bam
getStats ${locTMP}combined.polished.bam > ${OPENdir}/HiC.bam.stats.txt &


#convert to sorted bed
bedtools bamtobed -i ${locTMP}combined.polished.bam > ${locTMP}combined.polished.bed
sort --parallel=$THREADS -k 4 ${locTMP}combined.polished.bed >${locTMP}combined.polished.sorted.bed

#run SALSA
if [[ -n $GFA ]]; then
  GFAcommand="--gfa $GFA"
else
  GFAcommand=""
fi

SALSA -a $assemblyFASTA -l ${assemblyFASTA}.fai $GFAcommand -b ${locTMP}combined.polished.sorted.bed -e GATC -o ${OPENdir}scaffolded-salsa-noSplit-new/ --clean no -p yes 


wait
exit
