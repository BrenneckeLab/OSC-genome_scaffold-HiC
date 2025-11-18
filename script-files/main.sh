#!/usr/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --time=24:00:00
#SBATCH --qos=medium
#SBATCH --mem=10g


hostname

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

printf "start of pipeline\n\n" >${OPENdir}time-log.txt

###################################################################################################
#setup-phase

TMPdirRAW=$TMPdir

#prepare assembly file
if [[ $GAPFILLonly != Y ]]; then
  seqkit fx2tab $assemblyFASTA | grep -v N | seqkit tab2fx >${TMPdir}assembly.noN.fa
  assemblyFASTA=${TMPdir}assembly.noN.fa
fi
#---------------------------------------------------------------------------------------------------------
#fix contig names
if [[ ! -s ${TMPdir}genome.nameFixed.fa ]]; then
seqkit fx2tab $assemblyFASTA |
  awk '{ 
    n=split($1,X,/\:/)
    print ">"X[1]"~"NR"\n"$NF 
  }' > ${TMPdir}genome.nameFixed.fa
fi
#!do not comment out
assemblyFASTA=${TMPdir}genome.nameFixed.fa
VARI=${VARI},assemblyFASTA=${assemblyFASTA}
#!do not comment out

###################################################################################################
#start read-preparation for TGS-GapCloser
#extract runs to analyze
RUNs=$(echo $subRUNs | tr '~' '\t')
nRUNs=$(echo $RUNs | tr ' ' '\t' | wc -w)
VARI="${VARI},nRUNs=${nRUNs}"
SLURMid=

###################################################################################################
#prepare nanopore reads for gap-filling
##is run in the background while HiC scaffolding is happening

#remove existing fastq-file and regenerate
if [[ -z $NANOPOREreads ]]; then
  if [[ ! -s ${TMPdir}combined.fasta || $FORCE == Y ]]; then
    rm -rf ${TMPdir}reads.fq.gz
    for RUN in $RUNs; do
      cat ${RAW}${RUN}/reads/${BCversion}/trimmed/* >>${TMPdir}reads.fq.gz
    done

    cd $LOG
    #filter Nanopore reads with a bias on length
    SLURMid=$(sbatch --parsable -o "%x.o.%A.txt" -e "%x.e.%A.txt" --job-name=filtlong --mem=50g --cpus-per-task=2 --wrap="SINGULARITYdir=$SINGULARITYdir; source ${SCRIPTdir}tools; filtlong --min_length 20000 --target_bases 9000000000 --length_weight 10 ${TMPdir}reads.fq.gz | seqkit fq2fa --line-width 0 >${TMPdir}combined.fasta ") &

    #plot statistics of filtered Nanopore reads
    sbatch --parsable -o "%x.o.%A.txt" -e "%x.e.%A.txt" --dependency=$SLURMid --job-name=nanoplot --mem=50g --cpus-per-task=4 --wrap="SINGULARITYdir=$SINGULARITYdir; source ${SCRIPTdir}tools; NanoPlot -t \$SLURM_CPUS_PER_TASK --loglength --readtype 1D -o ${OPENdir}NanoPlot --format png --fasta ${TMPdir}combined.fasta  "
  fi

  NANOPOREreads=${TMPdir}combined.fasta
fi
VARI=${VARI},NANOPOREreads=$NANOPOREreads

#prepare raw-HiC data
locTMP=${TMPdir}prepareHiC/
mkdir -p $locTMP
echo $locTMP

if [[ ! -s ${locTMP}raw.1.fa || $FORCE == Y ]]; then
  #split paired end data into individual files
  #!remove head
  seqkit fx2tab $rawFASTA |
  mawk -v TMP=$locTMP '{
      sub("^>","",$1)
      if($1~"/1"){
        print ">"$1"\n"$2 > TMP "raw.1.fa"
      }else{
        print ">"$1"\n"$2 > TMP "raw.2.fa"
      }
  }'
fi
wait

###################################################################################################
###################################################################################################
if [[ ! -s ${OPENdir}correctedSALSA.agp ]]; then
  #!hacked in the manually split genome 
  #!still need to transfer the code to generate it also int this version
  #change to misassembly split genome
  assemblyFASTA=${TMPdir}genome.split.fa
  VARI=$VARI,assemblyFASTA=$assemblyFASTA

  ###################################################################################################
  #split misassemblies
  if [[ ! -s ${TMPdir}genome.split.fa ]]; then
    #detect misassemblies
    COMMAND="${SCRIPTdir}detectMisassemblies.sh"
    
    if [[ $COMPUTING == C ]]; then
      sbatch --dependency=$SLURMid --wait $COMMAND ${VARI}
    else
      #Simulate arrayID by adding it manually to VARI
      VARI="${VARI}"
      $COMMAND ${VARI}
    fi
  fi
  
  #change assembly to split genome
  assemblyFASTA=${TMPdir}genome.split.fa
  VARI=$VARI,assemblyFASTA=$assemblyFASTA
  

  #---------------------------------------------------------------------------------------------------------
  #create bwa index if it does not exist
  if [[ ! -f ${locTMP}bwa-index.amb ]]; then
    bwa index -p ${locTMP}bwa-index ${assemblyFASTA}
  fi

  #map HiC reads using bwa
  if [[ ! -s ${locTMP}mate_R1.bam ]]; then
      sbatch --wait --parsable --array=1-2 -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --job-name=map-bwa --mem=65g --cpus-per-task=38 --wrap="
    SINGULARITYdir=$SINGULARITYdir; THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 - 4 ));
    source ${SCRIPTdir}tools
    bwa mem -A1 -B4  -E50 -L0 -t \$THREADS ${locTMP}bwa-index ${locTMP}raw.\${SLURM_ARRAY_TASK_ID}.fa | samtools view -@ 4 -Shb - > ${locTMP}mate_R\${SLURM_ARRAY_TASK_ID}.bam"

  fi 

  #---------------------------------------------------------------------------------------------------------
  #run HiC analysis and scaffolding

  CUTOFF=SALSA
  BIN=detectMisassembly

  VARI=${VARI},CUTOFF=${CUTOFF},BIN=${BIN}

  COMMAND="${SCRIPTdir}process_HiC.sh"

  if [[ $COMPUTING == C ]]; then
    sbatch --wait $COMMAND ${VARI}
  else
    #Simulate arrayID by adding it manually to VARI
    VARI="${VARI}"
    $COMMAND ${VARI}
  fi

  cp ${OPENdir}scaffolded-salsa-noSplit-new/scaffolds_FINAL.fasta ${OPENdir}scaffolded-salsa-noSplit.fa
  FILLfile=${OPENdir}scaffolded-salsa-noSplit.fa

  COMMAND="annotate_assembly.sh -A $FILLfile -N ${assemblyVERSION}_scaffolded-${CUTOFF}-noSplit -r ${NANOPOREreads} -DHF"
  echo $COMMAND
  eval $COMMAND
  printf "\n\n\n please review the scaffolding \n in case there are problems please safe the modified agp file under exactly this path and rerun the pipeline:
  ${OPENdir}correctedSALSA.agp\n\n"
else
  agpcheck ${OPENdir}correctedSALSA.agp
  agp2fasta ${OPENdir}correctedSALSA.agp  ${TMPdir}genome.nameFixed.fa > ${OPENdir}SALSA_final.corrected.fa

  COMMAND="annotate_assembly.sh -A ${OPENdir}SALSA_final.corrected.fa -N ${assemblyVERSION}_scaffolded-SALSA-noSplit-corrected -r ${NANOPOREreads} -DHF"
  echo $COMMAND
  eval $COMMAND

fi

exit


sbatch --wait --parsable -o "%x.o.%A.txt" -e "%x.e.%A.txt" --job-name=map-convert_to_hic --mem=65g --cpus-per-task=4 --wrap="
  SINGULARITYdir=$SINGULARITYdir; THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2  ));
  source ${SCRIPTdir}tools
  convert_to_hic  ${OPENdir}scaffolded-salsa-noSplit-new/
"
exit


