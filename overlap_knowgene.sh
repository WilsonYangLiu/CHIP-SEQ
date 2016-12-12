# get kgID, gene symbol pairs
#grep -v "#" hg19.kgXref | awk -vOFS='\t' -vFS='\t' '{print $1,$5}' > kg-symbol

# define confident (FDR<=0.01%)
FDR=1
for sample in H3K56ac H3K56HMG; do
  # confident peaks
  grep -v "#" ${sample}_peaks.xls | awk -vOFS='\t' -vFDR=$FDR '(NR>2 && $9<=FDR){if($2<1){$2=1};print $1,$2,$3,$5,$7,$8,$9}' > ${sample}_confident
done

# get genes that each sample overlaped
for sample in H3K56ac H3K56HMG; do
  # option: -u (unique) Report the mere presence of any overlapping features in A
  awk -vOFS='\t' '{$2=$2+$4;$3=$2+1;print $0}' ${sample}_confident | intersectBed -a ./hg19_gene -b stdin -u | awk -vOFS='\t' '{print $4}' | sort -u > ${sample}_gene
done


