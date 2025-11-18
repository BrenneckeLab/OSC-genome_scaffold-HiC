#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=20g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=5:00:00

hostname
set -u

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase

#create path variables
locTMP=${TMPdir}fill_gaps/

#if rerun is requested wipe content of locTMP
if [[ $RERUN == Y ]]; then
  rm -rf ${locTMP}*
fi

#create directories
mkdir $locTMP

#load tools
source ${SCRIPTdir}tools

#calculate THREADS
THREADS=$(($SLURM_CPUS_PER_TASK * 2))
echo $THREADS
###################################################################################################
#prepare analysis

#!set assembly variable - even if file only copied in next step
currASSEMBLY=${locTMP}scaffolds.fa

#only run if fresh gap-filling run or if a re-run is requested
#!took away if statement as it is problematic if different scaffolding settings are tested
#!better to always run longer script here than to risk looking at the wrong version
if [[ ! -s ${locTMP}all-gaps.ID.seq.txt ]]; then
  #copy genome to local TMP
  rm -rf ${locTMP}scaffolds.fa*
  rm -rf ${OPENdir}gap-filling.log.txt

  #convert to single-line fasta and remove extra name-tag entires
  seqkit seq --only-id --upper-case --line-width 0 $FILLfile | seqkit rename --line-width 0 >${locTMP}scaffolds.fa
  cp $FILLfile ${locTMP}/
  #convert genome to 2bit
  faToTwoBit $currASSEMBLY ${locTMP}scaffolds.2bit

  #create choromosome size file
  seqkit fx2tab -n -l $currASSEMBLY | sort -k1,1 | tr -s '\t' '\t' >${locTMP}chrom.sizes

  #created file containing all gaps with names and identification-sequence
  twoBitInfo ${locTMP}scaffolds.2bit -nBed stdout | grep scaffold >${locTMP}all-gaps.bed

  #extend gap-locations by 1kb
  awk -v OFS="\t" '{
    $2=$2-1000
    $3=$3+1000
    $4=$1"::"NR
    print 
  }' ${locTMP}all-gaps.bed >${locTMP}all-gaps.extended.bed

  #extract gap-sequences for later gap-identification
  bedtools getfasta -tab -name -fi $currASSEMBLY -fo - -bed ${locTMP}all-gaps.extended.bed >${locTMP}all-gaps.ID.seq.txt
fi

#---------------------------------------------------------------------------------------------------------
    #create coordinates for flanking regions
    awk -v OFS="\t" -v SIZE=15000 -v CHROMsizes=${locTMP}chrom.sizes -v locTMP=$locTMP  '
    BEGIN{
      while((getline LINE < CHROMsizes) > 0) {
        split(LINE,splitLINE,/\t| /)
        CHRsize[splitLINE[1]]=splitLINE[2]
      }
    }
    {
      if($2<SIZE){
        print $1,0,$2,$1"::"NR":!:US"
      }else{
        print $1,$2-SIZE,$2,$1"::"NR":!:US"
      }
      if($3+SIZE > CHRsize[$1]){
        print $1,$3,CHRsize[$1],$1"::"NR":!:DS"
      }else{
        print $1,$3,$3+SIZE,$1"::"NR":!:DS"
      }
    }' ${locTMP}all-gaps.bed >${locTMP}flanking-regions.bed

    #extract sequence for flanking regions
    bedtools getfasta -fo ${locTMP}flanking-regions.fa -fi $currASSEMBLY -bed ${locTMP}flanking-regions.bed -name

    #map reads onto flanking regions
    minimap2 -2 -Q -x map-ont --secondary=no -t 60 ${locTMP}flanking-regions.fa ${NANOPOREreads} >${locTMP}reads_mapped_to_flankingRegions.paf

#---------------------------------------------------------------------------------------------------------
#!get gap-coordinate in current genome version
twoBitInfo ${locTMP}scaffolds.2bit -nBed stdout

#get the first gap
currGAP=$(twoBitInfo ${locTMP}scaffolds.2bit -nBed stdout | grep scaffold | head -n1 | tail -n 1)
currGAP=$(grep "18769095" ${locTMP}all-gaps.bed )
echo $currGAP

#extend gap-locations by 1kb
awk -v OFS="\t" '{
  $2=$2-1000
  $3=$3+1000
  $4=$1"::"NR
  print 
}' <(echo $currGAP) >${locTMP}currGAP.extended.bed

#extract gap-sequences for later gap-identification
bedtools getfasta -tab -name -fi $currASSEMBLY -fo - -bed ${locTMP}currGAP.extended.bed >${locTMP}currGAP.ID.seq.txt

#extract GAP ID based on sequence
currGAPseq=$(cut -f 2 ${locTMP}currGAP.ID.seq.txt)
GAPid=$(grep $currGAPseq ${locTMP}all-gaps.ID.seq.txt | cut -f 1)

#---------------------------------------------------------------------------------------------------------
#prepare loop
#clear
rm -rf ${locTMP}replaced.all.fa

#initiate skip-ahead variable
SKIP=0
allGAPprocessed=1
splitCONTIG=

###################################################################################################
#loop through all the gaps
while [[ ! -z $currGAP ]]; do

  printf "\n\n\nprocessing $GAPid \n\n\n"
  echo $GAPid >>${OPENdir}gap-filling.log.txt
  echo $currGAP >>${OPENdir}gap-filling.log.txt

  CHR=$(echo $currGAP | tr ' ' '\t' | cut -f 1)
  START=$(echo $currGAP | tr ' ' '\t' | cut -f 2)
  STOP=$(echo $currGAP | tr ' ' '\t' | cut -f 3)

  currTMP=${locTMP}${GAPid}/
  rm -rf ${currTMP}
  mkdir $currTMP

  #only run identifiaction of flanking regions and associated reads if not already existing
  if [[ ! -s ${currTMP}/reads.fa ]]; then
    #create coordinates for flanking regions
    grep $GAPid ${locTMP}flanking-regions.bed > ${currTMP}flanking-regions.bed

    #extract sequence for flanking regions
    grep -A 1 $GAPid ${locTMP}flanking-regions.fa > ${currTMP}flanking-regions.fa

    #map reads onto flanking regions
    grep $GAPid ${locTMP}reads_mapped_to_flankingRegions.paf > ${currTMP}reads_mapped_to_flankingRegions.paf

    #extract IDs of reads overlapping the flanking regions
    awk -v OFS="\t" '
    {
        if($9-$8 >$7*0.8){
          print $1
        }
    }' ${currTMP}reads_mapped_to_flankingRegions.paf |
      sort | uniq >${currTMP}IDs.txt

    #extract gap-reads from the full fasta file
    seqkit grep --line-width 0 -n -f ${currTMP}/IDs.txt ${NANOPOREreads} >${currTMP}/reads.fa
  fi

  #determine if both fanking regions contained within a single read
  rm -rf ${currTMP}bridging-ID.longest.txt
  awk -v OFS="\t" '{
    if($9-$8 > $7*0.8 && $4+5000 < $2){
      split($6,splitNAME,/:!:/)
      X[$1]["COUNT"]+=1
      X[$1][splitNAME[2]]=$0
    }
  }
  END{
    for(READ in X){
      if("US" in X[READ] && "DS" in X[READ] && X[READ]["COUNT"]==2 ){
        n=split(READ,splitREAD,/::|=/)
        for(i=1; i<=n; i++0){if(splitREAD[i]=="length") READlength=splitREAD[i+1]}
        print READ,READlength
      }
    }
  }' ${currTMP}reads_mapped_to_flankingRegions.paf | sort -k2,2nr | head -n 1 | cut -f1 >${currTMP}bridging-ID.longest.txt

  #preset empty SCAFFstat variable
  SCAFFstat=

  #if a bridging read exists use this for polishing and then gapfilling
  if [[ -s ${currTMP}bridging-ID.longest.txt ]]; then
    echo bridging read >>${OPENdir}gap-filling.log.txt
    

    #get sequence for longest bridging read
    seqkit grep --line-width 0 -n -f ${currTMP}bridging-ID.longest.txt ${NANOPOREreads} >${currTMP}/longest-read.fa
    seqkit grep --line-width 0 -v -n -f ${currTMP}bridging-ID.longest.txt ${NANOPOREreads} >${currTMP}/reads.exclLongest.fa
    nREADS=$(cat ${currTMP}/longest-read.fa | wc -l)

    if [[ $nREADS -le 2 ]]; then

      #map all reads to the bridging read
      minimap2 -2 -Q -x map-ont --secondary=no -t $THREADS ${currTMP}/longest-read.fa ${currTMP}/reads.exclLongest.fa >${currTMP}reads_mapped_bridging.paf

      #polish bridging read
      racon -m 8 -x -6 -g -8 -w 500 -t $THREADS ${currTMP}/reads.exclLongest.fa ${currTMP}reads_mapped_bridging.paf ${currTMP}/longest-read.fa >${currTMP}bridging.polished.fa

      #generate contig.fa required for replacing the gap in the genome
      rm -rf ${currTMP}contig.fa
      cp ${currTMP}bridging.polished.fa ${currTMP}contig.fa
    fi

  fi

  if [[ ! -s ${currTMP}contig.fa ]]; then
    echo assemble reads >>${OPENdir}gap-filling.log.txt

    #prepare reads for assembly
    cat ${currTMP}/reads.fa | gzip >${currTMP}/reads.fa.gz

    #assemble reads
    wtdbg2 -x ont -g 100k -i ${currTMP}/reads.fa.gz -t 16 -fo ${currTMP}wtdbg2-out
    wtpoa -t 16 -i ${currTMP}wtdbg2-out.ctg.lay.gz -fo ${currTMP}dbg.raw.fa
    minimap2 -2 -Q -x map-ont --secondary=no -t $THREADS ${currTMP}dbg.raw.fa ${currTMP}/reads.fa >${currTMP}reads_mapped_wtcbg-contig.paf

    mkdir ${currTMP}out-corr/
    racon -m 8 -x -6 -g -8 -w 500 -t $THREADS ${currTMP}/reads.fa ${currTMP}reads_mapped_wtcbg-contig.paf ${currTMP}dbg.raw.fa >${currTMP}out-corr/assembly.fasta

    #test if both flanking regions contained in a single contig
    minimap2 -2 -Q -x map-ont --secondary=no -t 60 ${currTMP}out-corr/assembly.fasta ${currTMP}flanking-regions.fa >${currTMP}flanking_mapped_to_contigs.paf

    #reset variable and splitCONTIG file for next round
    OKcontig=
    rm -rf ${currTMP}splitCONTIG.bed
    #test if a contig or scaffold containing both flanking regions exist
    OKcontig=$(awk -v OFS="\t" -v currTMP=$currTMP '
    {  
        if($4-$3 >$2*0.8){
          split($1,splitNAME,/:!:/)
          X[$6][splitNAME[2]]=$0
        }
    }
    END{
      for(i in X){
        if("US" in X[i] && "DS" in X[i]){
          OK=i
          nOK+=1
        }
      }
      if(nOK == 1){
        print OK
      }else{
        for(i in X){
          if("US" in X[i] ){
            nUS+=1
            US=i
            split(X[i]["US"], splitLINE,/ |\t/)
            USstrand=splitLINE[5]
            USstart=splitLINE[8]
            USlength=splitLINE[7]
          }
          if("DS" in X[i] ){
            nDS+=1
            DS=i
            split(X[i]["DS"], splitLINE,/ |\t/)
            DSstrand=splitLINE[5]
            DSstart=splitLINE[8]
            DSlength=splitLINE[7]
          }
        }
        if(nUS == 1 && nDS == 1 ){
          if(USstart < DSstart){
            print US,0,USlength,"UScontig",0,USstrand > currTMP "splitCONTIG.bed"
            print DS,0,DSlength,"DScontig",0,DSstrand > currTMP "splitCONTIG.bed"
          }else{
            print DS,0,DSlength,"DScontig",0,DSstrand > currTMP "splitCONTIG.bed"
            print US,0,USlength,"UScontig",0,USstrand > currTMP "splitCONTIG.bed"
          }
        }else{
          print "not exactly 1 contig or scaffold containing both flanking regions or clean split-contigs" > currTMP "evaluation.log"
        }
      }
    }' ${currTMP}flanking_mapped_to_contigs.paf)

    rm -rf ${currTMP}contig.fa
    if [[ ! -z $OKcontig ]]; then
      echo $OKcontig contains both flanking regions >>${OPENdir}gap-filling.log.txt
      seqkit grep --line-width 0 -p $OKcontig ${currTMP}/out-corr/assembly.fasta >${currTMP}contig.fa

      #set variable if scaffold used for gapfilling
      SCAFFstat=$(grep -v ">" ${currTMP}contig.fa | grep N)
    elif [[ -s ${currTMP}splitCONTIG.bed ]]; then
      splitCONTIG=Y
      echo split contig available - fusing with Ns >>${OPENdir}gap-filling.log.txt
      bedtools getfasta -tab -s -name -fi ${currTMP}/out-corr/assembly.fasta -fo ${currTMP}splitCONTIG.tab -bed ${currTMP}splitCONTIG.bed
      awk -v OFS="\t" '
      {
        if(NR==1){
          print ">splitCONTIG"
          X=$NF"NNNNNNNNNNNNNNN"
        }else{
          X=X$NF
        }
      }
      END{
        print X
      }' ${currTMP}splitCONTIG.tab >${currTMP}contig.fa
    fi

  fi

  #!if only single contig
  if [[ -s ${currTMP}contig.fa ]]; then
    #verify that both flanking regions are contained in the contig
    minimap2 -2 -Q -x map-ont -t 60 ${currTMP}contig.fa ${currTMP}flanking-regions.fa >${currTMP}flanking_to_contig.paf

    #test output if only 2 valid flanking region mappings are present
    awk -v OUTFILE=${OPENdir}gap-filling.log.txt '
    {
        if($4-$3 > $2*0.8 ){
          split($1,splitNAME,/:!:/)
          X["COUNT"]+=1
          X["FLANK"][splitNAME[2]]=$0 }
      }
      END{
        if(X["COUNT"]==2) {
          if("US" in X["FLANK"] && "DS" in X["FLANK"]){
            print "OK"
          }else{
            print "not both flanks present" >> OUTFILE
          }
        }else{
          if(X["COUNT"]<2){
            print "less than 2 flanking region-mappings" >> OUTFILE
          }else{
            print "more than 2 flanking region-mappings" >> OUTFILE
          }
        }
      }' ${currTMP}flanking_to_contig.paf >${currTMP}evaluation.log

    #!evaluate status
    STATUS=$(cat ${currTMP}evaluation.log)
    echo $STATUS

    #!STATUS test
    if [[ $STATUS == OK ]]; then
      #delete error.txt
      rm -rf ${currTMP}error.txt

      #map contig to the genome
      minimap2 -2 -Q -x map-ont --secondary=no -t 60 $currASSEMBLY ${currTMP}contig.fa >${currTMP}contig_to_scaffold.paf
      #check if everything is ok and if so output the bed-files required for sequence replacement
      awk -v OFS="\t" -v SIZE=15000 -v currTMP=$currTMP -v CHR=$CHR -v START=$START -v STOP=$STOP '{
          if($4-$3>SIZE && $6 == CHR ){
            X+=1
            Y[NR]=$0  
          }
        }
        END{
          #only alow up to 2 mappings (either spanning gap or split alignment)
          if(X <= 2){
            contigSTART=10000000000
            contigSTOP=0
            chromSTART=1000000000000000
            chromSTOP=0
            for( i in Y){
              split(Y[i],splitLINE,/ |\t/)
              if(i == 1){
                STRAND=splitLINE[5]
              }
              if( STRAND == splitLINE[5]){
                CONTIG=splitLINE[1]
                if(splitLINE[3]<contigSTART){contigSTART=splitLINE[3]}
                if(splitLINE[4]>contigSTOP){contigSTOP=splitLINE[4]}
                CHROM=splitLINE[6]
                print Y[i]
                if(splitLINE[8]<chromSTART){chromSTART=splitLINE[8]}
                if(splitLINE[9]>chromSTOP){chromSTOP=splitLINE[9]}

              }else{
                print "different strands" > currTMP "error.txt"
                for(i in Y){
                  print Y[i]> currTMP "error.txt"
                }
                exit 1
              }
            }
            if( CHROM==CHR && chromSTART<START && chromSTOP>STOP){
              print CONTIG,contigSTART,contigSTOP,"contig",1,STRAND > currTMP "contig.bed"
              print CHROM,chromSTART,chromSTOP,"chrom",1,"+" > currTMP "chrom.bed"
            }else{
              print "ohoh",CHROM"=="CHR,chromSTART"<"START,chromSTOP">"STOP > currTMP "error.txt"
            }

          }else{
            print  "not 2 mappings" > currTMP "error.txt"
            for(i in Y){
              print Y[i]> currTMP "error.txt"
            }
            exit 1
          }
        } ' ${currTMP}contig_to_scaffold.paf

      #if everything is fine replace the sequence
      if [[ ! -s ${currTMP}error.txt ]]; then
        #extract contig sequence
        bedtools getfasta -s -tab -fi ${currTMP}contig.fa -fo ${currTMP}replacement.tab -bed ${currTMP}contig.bed
        bedtools getfasta -s -tab -fi $currASSEMBLY -fo ${currTMP}toReplace.tab -bed ${currTMP}chrom.bed

        toREPLACE=$(cat ${currTMP}toReplace.tab | cut -f 2)
        cat ${currTMP}replacement.tab | cut -f 2 >${currTMP}contig.seq

        seqkit fx2tab --line-width 0 $currASSEMBLY |
          awk -v contigSEQ=${currTMP}contig.seq -v chromBED=${currTMP}chrom.bed '
          BEGIN{
            while((getline LINE < chromBED) > 0) {
              split(LINE,splitLINE,/ |\t/)
              CHROM=splitLINE[1]
              START=splitLINE[2]
              STOP=splitLINE[3]
            }
            while((getline LINE < contigSEQ) > 0) {
              SEQ=LINE
            }
          }
          {
            if($1 == CHROM){
              print ">"$1"\n" substr($NF,1,START) SEQ substr($NF,STOP,length($NF))
              #print length(SEQ)
            }else{
              print ">"$1 "\n" $NF 
              #a=b
            }
          }' >${currTMP}genome.fa

        rm -rf ${currASSEMBLY}
        rm -rf ${currASSEMBLY}.fai
        cp ${currTMP}genome.fa $currASSEMBLY

        #add inserted sequence to file for later mapping to get coordinates
        echo ">"$currGAP >>${locTMP}replaced.all.fa
        cut -f 2 ${currTMP}replacement.tab >>${locTMP}replaced.all.fa

        printf "succesfully replaced gap \n\n\n" >>${OPENdir}gap-filling.log.txt

        if [[ $splitCONTIG == Y ]]; then
          SKIP=$(($SKIP + 1))
        fi
        splitCONTIG=
      else
        cat ${currTMP}error.txt >>${OPENdir}gap-filling.log.txt
        echo $currGAP - error.txt present
        SKIP=$(($SKIP + 1))
        printf "contig does not match genome position \n\n\n" >>${OPENdir}gap-filling.log.txt
      fi
    else
      echo $currGAP - OK status not validated
      SKIP=$(($SKIP + 1))
      printf "contig status not OK \n\n\n" >>${OPENdir}gap-filling.log.txt
    fi

  else
    echo $currGAP - multiple contigs
    SKIP=$(($SKIP + 1))
    printf "multiple contigs \n\n\n" >>${OPENdir}gap-filling.log.txt
  fi

  #convert genome to 2bit
  rm -rf ${locTMP}scaffolds.2bit
  faToTwoBit $currASSEMBLY ${locTMP}scaffolds.2bit

  #create choromosome size file
  seqkit fx2tab -n -l $currASSEMBLY | sort -k1,1 | tr -s '\t' '\t' >${locTMP}chrom.sizes

  #if scaffold was used for filling add new flanking regions to ID-file
  if [[ ! -z $SCAFFstat ]]; then
    #created file containing all gaps with names and identification-sequence
    twoBitInfo ${locTMP}scaffolds.2bit -nBed stdout | grep scaffold | head -n 1 >${locTMP}scaff-gap.bed

    #extend gap-locations by 1kb
    awk -v OFS="\t" -v GAPid=$GAPid '{
    $2=$2-1000
    $3=$3+1000
    $4=GAPid"::"NR"::sacffold"
    print 
  }' ${locTMP}scaff-gap.bed >${locTMP}scaff-gap.extended.bed

    #extract gap-sequences for later gap-identification
    bedtools getfasta -tab -name -fi $currASSEMBLY -fo - -bed scaff-gap.extended.bed >>${locTMP}all-gaps.ID.seq.txt
  fi

  #!get next gap location
  #add skip ahead to extract the next gap to process
  nGAP=$((1 + $SKIP))
  allGAPprocessed=$(( $allGAPprocessed + 1 ))
  nGAPtotal=$(cat ${locTMP}all-gaps.ID.seq.txt | wc -l)

  currGAP=

  if [[ $allGAPprocessed -le $nGAPtotal ]]; then
    currGAP=$(twoBitInfo ${locTMP}scaffolds.2bit -nBed stdout | grep scaffold | head -n $nGAP | tail -n 1)
    if [[ ! -z $currGAP ]]; then

      #extend gap-locations by 1kb
      awk -v OFS="\t" '{
      $2=$2-1000
      $3=$3+1000
      $4=$1"::"NR
      print 
    }' <(echo $currGAP) >${locTMP}currGAP.extended.bed

      #extract gap-sequences for later gap-identification
      bedtools getfasta -tab -name -fi $currASSEMBLY -fo - -bed ${locTMP}currGAP.extended.bed >${locTMP}currGAP.ID.seq.txt

      #extract GAP ID based on sequence
      currGAPseq=$(cut -f 2 ${locTMP}currGAP.ID.seq.txt)
      GAPid=$(grep $currGAPseq ${locTMP}all-gaps.ID.seq.txt | cut -f 1)
    fi
  fi
done

cp $currASSEMBLY ${OPENdir}gap-filled-${CUTOFF}-${BIN}.fa
exit
