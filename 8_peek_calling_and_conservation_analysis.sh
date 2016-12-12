# peek calling. Using MACS call peeks for each IP sample and its corresponding input control sample(cell_lysate)
:<<BLOCK
for pair in H3K56ac-cell_lysate H3K56HMG-cell_lysate IgG-cell_lysate; do
  echo -en $pair"\t"
  chip=$(echo $pair | sed 's/-.*//')
  input=$(echo $pair | sed 's/.*-//')
  # Run MACS
  GEN_SIZE=hs
  PVALUE=1e-5
  SHIFTSIZE=73
  
  macs14 -t ${chip}.bam -c ${input}.bam --pvalue=$PVALUE --gsize=$GEN_SIZE --format=BAM --name=${chip} --nomodel --shiftsize=$SHIFTSIZE 2> ${pair}_MACS.log

  # check warnings
  grep "WARNING" ${pair}_MACS.log
  # rm intermediate files
  
done
BLOCK

# number of peaks at different FDR threadholds
echo -e "FDR\tall\t5\t1\t0"
for sample in H3K56ac H3K56HMG IgG; do
  echo -en ${sample}-cl
  for fdr in 100 5 1 0; do
    echo -en "\t"$(grep -v "#" ${sample}_peaks.xls | awk -v FDR=$fdr 'NR>2 && $9<=FDR' | wc -l)
  done
  echo	# newline
done 

# define confident (FDR), enriched regions (p >= 1e-5) and control peaks
FDR=1
for sample in H3K56ac H3K56HMG IgG; do
  # confident peaks
  grep -v "#" ${sample}_peaks.xls | awk -vOFS='\t' -vFDR=$FDR '(NR>2 && $9<=FDR){if($2<1){$2=1};print $1,$2,$3,$5,$7,$8,$9}' > ${sample}_confident
  # Regions with significant enrichment
  grep -v "#" ${sample}_peaks.xls | awk -vOFS='\t' '(NR>2){if($2<1){$2=1};print $1,$2,$3,$5,$7,$8,$9}' > ${sample}_enrichment
  # Control peaks
  shuffleBed -i ${sample}_enrichment -g hg19.genome -chrom | sort -k1,1 -k2,2n > ${sample}_control
done

# Peak visualization
for sample in H3K56ac H3K56HMG IgG; do
  # Create BED files
  (echo -e "track name=\"${sample}_confident_peaks\" description=\"${sample}_confident_peaks\" visibility=2"
  sort -k5,5gr ${sample}_confident | awk -vOFS='\t' '{print $1,$2,$3,"PEAK_"NR,$5,"."}' | sort -k1,1 -k2,2n) | gzip > ${sample}_confident.bed.gz
  (echo -e "track name=\"${sample}_enrichment_peaks\" description=\"${sample}_enrichment_peaks\" visibility=2"
  sort -k5,5gr ${sample}_enrichment | awk -vOFS='\t' '{print $1,$2,$3,"PEAK_"NR,$5,"."}' | sort -k1,1 -k2,2n) | gzip > ${sample}_enrichment.bed.gz
done

# Caculate a conservation rate between reference and sample. Also caculate the conservation rate for control peaks.
# Note that if the number of peaks in very different between two conditions, the rate of binding conservation depend on which condition is chosen as the reference sample.
# Overlap summit of reference confident peaks with sample enriched regions and reference control peaks
for pair in H3K56ac-H3K56HMG H3K56HMG-H3K56ac; do
  reference=$(echo $pair | sed 's/-.*//')
  sample=$(echo $pair | sed 's/.*-//')
  TOTAL=$(cat ${reference}_confident | wc -l)
  echo -e "Condition(ref-sample)\t# of total\t# of overlap\t% of overlap"
  awk -vOFS='\t' '{$2=$2+$4;$3=$2+1;print $0}' ${reference}_confident | intersectBed -a ${sample}_enrichment -b stdin  | wc -l | awk -vOFS='\t' -vref=$reference -vsample=$sample -vTOTAL=$TOTAL '{print ref"-"sample,TOTAL,$1,$1*100/TOTAL}'
  awk -vOFS='\t' '{$2=$2+$4;$3=$2+1;print $0}' ${reference}_confident | intersectBed -a ${reference}_control -b stdin | wc -l | awk -vOFS='\t' -vref=$reference -vTOTAL=$TOTAL '{print ref"-ref_control",TOTAL,$1,$1*100/TOTAL}'
done


