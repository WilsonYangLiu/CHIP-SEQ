#####
## Example code
## Bardet AF, He Q, Zeitlinger J. & Stark A
## A computational pipeline for comparative ChIP-seq analyses
## Nat Protoc. 2011 Dec 15;7(1):45-61.
#####

## Test dataset can be downloaded from the data section of www.starklab.org
# All commands are BASH commands

## Step 1: Sequence quality check
for sample in chip_dmel input_dmel chip_dyak input_dyak; do
    # FASTX Statistics
    fastx_quality_stats -i <(gunzip -c ${sample}.fastq.gz) -o ${sample}_stats.txt
    # FASTX quality score
    fastq_quality_boxplot_graph.sh -i ${sample}_stats.txt -o ${sample}_quality.png -t ${sample}
    # FASTX nucleotide distribution
    fastx_nucleotide_distribution_graph.sh -i ${sample}_stats.txt -o ${sample}_nuc.png -t ${sample}
    # Remove intermediate file
    rm ${sample}_stats.txt
done

## Step 2: Raw reads count
for sample in chip_dmel input_dmel chip_dyak input_dyak; do
    echo -en $sample"\t"
    # Number of unique reads and most repeated read
    gunzip -c ${sample}.fastq.gz | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max){max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'
done

## Step 3: Read length check
for sample in chip_dmel input_dmel chip_dyak input_dyak; do
    echo -en $sample"\t"
    # Read length
    gunzip -c ${sample}.fastq.gz | awk '((NR-2)%4==0){count[length($1)]++}END{for(len in count){print len}}'
    # Truncate longer reads to 36bp (if necessary)
    LEN=36
    gunzip -c ${sample}.fastq.gz | awk -vLEN=$LEN '{if((NR-2)%2==0){print substr($1,1,LEN)}else{print $0}}' | gzip > ${sample}_36bp.fastq.gz
done

## Step 4: Read mapping
# Map reads
for sample in chip_dmel input_dmel; do
    gunzip -c ${sample}_36bp.fastq.gz | bowtie -q -m 1 -v 3 --sam --best --strata bowtie_index_dm3/dm3 - > ${sample}.sam
done
for sample in chip_dyak input_dyak; do
    gunzip -c ${sample}_36bp.fastq.gz | bowtie -q -m 1 -v 3 --sam --best --strata bowtie_index_droYak2/droYak2 - > ${sample}.sam
done
for sample in chip_dmel input_dmel chip_dyak input_dyak; do
    # Convert file from SAM to BAM format
    samtools view -Sb ${sample}.sam > ${sample}_nonSorted.bam
    # Sort BAM file
    samtools sort ${sample}_nonSorted.bam ${sample}
    # Create index file (BAI)
    samtools index ${sample}.bam
    # Revove intermediate files
    rm ${sample}.sam ${sample}_nonSorted.bam
done

## Step 5: Mapped reads counts
for sample in chip_dmel input_dmel chip_dyak input_dyak; do
    echo -en $sample"\t"
    # Number of raw reads
    raw=$(samtools view ${sample}.bam | wc -l)
    # Number of raw, unique and most repeated reads
    bamToBed -i ${sample}.bam | awk -vRAW=$raw '{coordinates=$1":"$2"-"$3;total++;count[coordinates]++}END{for(coordinates in count){if(!max||count[coordinates]>max){max=count[coordinates];maxCoor= coordinates};if(count[coordinates]==1){unique++}};print RAW,total,total*100/RAW,unique,unique*100/total,maxCoor,count[maxCoor],count[maxCoor]*100/total}'
    # Total and top 10 of non-mapped reads
    samtools view -f 0x0004 ${sample}.bam | awk '{read=$10;total++;count[read]++}END{print "Total_non-mapped_reads",total;for(read in count){print read,count[read]+0}}' | sort -k2,2nr | head -11
done

## Step 6: Read density visualization
for sample in chip_dmel input_dmel; do
    EXTEND=150
    # Number of reads
    librarySize=$(samtools idxstats ${sample}.bam | awk '{total+=$3} END{print total}')
    # Create density file: extend reads, calculate read density at each position and normalize the library size to 1 million reads
    bamToBed -i ${sample}.bam | awk -vCHROM="dm3.chrom.sizes" -vEXTEND=$EXTEND -vOFS='\t' 'BEGIN{while(getline<CHROM){chromSize[$1]=$2}}{chrom=$1;start=$2;end=$3;strand=$6;if(strand=="+"){end=start+EXTEND;if(end>chromSize[chrom]){end=chromSize[chrom]}};if(strand=="-"){start=end-EXTEND;if(start<1){start=1}};print chrom,start,end}' | sort -k1,1 -k2,2n | genomeCoverageBed -i stdin -g dm3.chrom.sizes -d | awk -vOFS='\t' -vSIZE=$librarySize '{print $1,$2,$2+1,$3*1000000/SIZE}' | gzip > ${sample}.density.gz
     # Create WIG file
    gunzip -c ${sample}.density.gz | awk -vOFS='\t' '($4!=0){if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | gzip > ${sample}.wig.gz
    # Create BigWig file
    wigToBigWig ${sample}.wig.gz dm3.chrom.sizes ${sample}.bw
    # Remove intermediate file
    rm ${sample}.wig.gz
done

## Box 1 Step 1: Read translation

# A. Standard processing on single core machine (long running time)
for sample in chip_dyak input_dyak; do
    # Translate the coordinates from genome to reference genome and keep information of where the reads came from in the genome in the name column of the BED file
    liftOver <(bamToBed -i ${sample}.bam | awk -vOFS='\t' '{$4=$1":"$2":"$3;print $0}') droYak2Todm3.over.chain ${sample}_dm3_tmp.bed ${sample}_dm3_lost.bed
done

# B. Parallel processing on multi-core machines
for sample in chip_dyak input_dyak; do
    split=5
    # Split big input file in 5 smaller files and keep information of where the reads came from in the genome in the name column of the BED file
    bamToBed -i ${sample}.bam | awk -vOFS='\t' -vSPLIT=$split -vFILE=${sample} '{$4=$1":"$2":"$3;print $0>(FILE"_"(NR%SPLIT)+1".bed")}'
    # Translate the coordinates from genome to reference genome
    for i in `seq 1 1 $split`; do 
	liftOver ${sample}_${i}.bed droYak2Todm3.over.chain ${sample}_${i}_dm3_tmp.bed ${sample}_${i}_dm3_lost.bed &
    done
done
# Merge output files
for sample in chip_dyak input_dyak; do
    sort -k1,1 -k2,2n ${sample}_*_dm3_tmp.bed > ${sample}_dm3_tmp.bed
    sort -k1,1 -k2,2n ${sample}_*_dm3_lost.bed > ${sample}_dm3_lost.bed
    rm ${sample}_[0-9]*.bed
done

## Box 1 Step 2: Translated reads counts
for sample in chip_dyak_dm3 input_dyak_dm3; do
    PERCENT=10
    # Remove reads which length changed by more than 10%
    awk -vPERCENT=$PERCENT '{split($4,COOR,":");lengthBefore=COOR[3]-COOR[2];lengthAfter=$3-$2;if(lengthAfter>(lengthBefore*(100-PERCENT)/100)&&lengthAfter<(lengthBefore*(100+PERCENT)/100)){print $0}}' ${sample}_tmp.bed | grep -v "chrU"> ${sample}.bed
    # Count number of translated reads
    wc -l ${sample}_tmp.bed ${sample}.bed
    # Convert BED to BAM file
    bedToBam -i ${sample}.bed -g dm3.chrom.sizes > ${sample}_nonSorted.bam
    # Sort BAM file
    samtools sort ${sample}_nonSorted.bam ${sample}
    # Create index file (BAI)
    samtools index ${sample}.bam
    # Remove intermediate files
    rm ${sample}_nonSorted.bam ${sample}_lost.bed ${sample}_tmp.bed ${sample}.bed
done

## Box 1 Step 3: Read density visualization
for sample in chip_dyak_dm3 input_dyak_dm3; do
    EXTEND=150
    # Number of reads
    librarySize=$(samtools idxstats ${sample}.bam | awk '{total+=$3} END{print total}')
    # Create density file: extend reads, calculate read density at each position and normalize the library size to 1 million reads
    bamToBed -i ${sample}.bam | awk -vCHROM="dm3.chrom.sizes" -vEXTEND=$EXTEND -vOFS='\t' 'BEGIN{while(getline<CHROM){chromSize[$1]=$2}}{chrom=$1;start=$2;end=$3;strand=$6;if(strand=="+"){end=start+EXTEND;if(end>chromSize[chrom]){end=chromSize[chrom]}};if(strand=="-"){start=end-EXTEND;if(start<1){start=1}};print chrom,start,end}' | sort -k1,1 -k2,2n | genomeCoverageBed -i stdin -g dm3.chrom.sizes -d | awk -vOFS='\t' -vSIZE=$librarySize '{print $1,$2,$2+1,$3*1000000/SIZE}' | gzip > ${sample}.density.gz
     # Create WIG file
    gunzip -c ${sample}.density.gz | awk -vOFS='\t' '($4!=0){if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | gzip > ${sample}.wig.gz
    # Create BigWig file
    wigToBigWig ${sample}.wig.gz dm3.chrom.sizes ${sample}.bw
    # Remove intermediate file
    rm ${sample}.wig.gz
done

## Step 7: Pearson correlation
for pair in chip_dmel-input_dmel chip_dyak_dm3-input_dyak_dm3 chip_dmel-chip_dyak_dm3; do
    echo -en $sample"\t"
    chip=$(echo $pair | sed 's/-.*//')
    input=$(echo $pair | sed 's/.*-//')
    paste <(gunzip -c ${chip}.density.gz) <(gunzip -c ${input}.density.gz) | awk '{if($2!=$6){exit 1};if($4!=0||$8!=0){print $4,$8}}' | correlation.awk
done

## Step 8: Peak calling
for pair in chip_dmel-input_dmel chip_dyak_dm3-input_dyak_dm3; do
    echo -en $pair"\t"
    chip=$(echo $pair | sed 's/-.*//')
    input=$(echo $pair | sed 's/.*-//')
    # Run MACS
    GEN_SIZE=$(awk '{size+=$2}END{print size}' dm3.chrom.sizes)
    READ_LEN=36
    PVALUE=1e-5
    MFOLD=4 # Maximum possible
    macs -t ${chip}.bam -c ${input}.bam --name=${pair}_macs_p05 --format=BAM --gsize=$GEN_SIZE --tsize=$READ_LEN --pvalue=$PVALUE --mfold=$MFOLD 2> ${pair}_macs_p05.log
    # Print genomic fragment length
    grep "# d = " ${pair}_macs_p05_peaks.xls | awk '{print $4}'
    # Check warnings
    grep "WARNING" ${pair}_macs_p05.log
    # Remove intermediate files
    rm ${pair}_macs_p05{.log,_model.r,_negative_peaks.xls,_peaks.bed}
done
# Re-run MACS with a lower mfold parameter if WARNING "Fewer paired peaks (X) than 1000!"
# Number of peaks at different FDR thresholds
(echo -e "FDR\tAll\t5\t1\t0"
for pair in chip_dmel-input_dmel chip_dyak_dm3-input_dyak_dm3; do
    echo -en $pair
    for fdr in 100 5 1 0; do
        echo -en "\t"$(grep -v "#" ${pair}_macs_p05_peaks.xls | awk -vFDR=$fdr '(NR>1&&$9<=FDR)' | wc -l)
    done
    echo
done)
# Define confident peaks (FDR), enriched regions (p-value<=10e-5) and control peaks
FDR=1
for pair in chip_dmel-input_dmel chip_dyak_dm3-input_dyak_dm3; do
    # Confident peaks
    grep -v "#" ${pair}_macs_p05_peaks.xls | awk -vOFS='\t' -vFDR=$FDR '(NR>1&&$9<=FDR){if($2<1){$2=1};print $1,$2,$3,$5,$7,$8,$9}' > ${pair}_macs_confident.txt
    # Regions with significant enrichment
    grep -v "#" ${pair}_macs_p05_peaks.xls | awk -vOFS='\t' '(NR>1){if($2<1){$2=1};print $1,$2,$3,$5,$7,$8,$9}' > ${pair}_macs_enrichment.txt
    # Control peaks
    shuffleBed -i ${pair}_macs_enrichment.txt -g dm3.chrom.sizes -chrom | sort -k1,1 -k2,2n > ${pair}_macs_control.txt
done

## Step 9: Peak visualization
for pair in chip_dmel-input_dmel chip_dyak_dm3-input_dyak_dm3; do
    # Create BED files
    (echo -e "track name=\"${pair}_confident_peaks\" description=\"${pair}_confident_peaks\" visibility=2"
    sort -k5,5gr ${pair}_macs_confident.txt | awk -vOFS='\t' '{print $1,$2,$3,"PEAK_"NR,$5,"."}' | sort -k1,1 -k2,2n) | gzip > ${pair}_macs_confident.bed.gz
    (echo -e "track name=\"${pair}_enriched_regions\" description=\"${pair}_enriched_regions\" visibility=2"
    sort -k5,5gr ${pair}_macs_enrichment.txt | awk -vOFS='\t' '{print $1,$2,$3,"PEAK_"NR,$5,"."}' | sort -k1,1 -k2,2n) | gzip > ${pair}_macs_enrichment.bed.gz
done

## Step 10: Peak conservation
reference=chip_dmel-input_dmel
sample=chip_dyak_dm3-input_dyak_dm3
# Overlap summit of reference confident peaks with sample enriched regions and reference control peaks
TOTAL=$(cat ${reference}_macs_confident.txt | wc -l)
awk -vOFS='\t' '{$2=$2+$4;$3=$2+1;print $0}' ${reference}_macs_confident.txt | intersectBed -a stdin -b ${sample}_macs_enrichment.txt | wc -l | awk -vTOTAL=$TOTAL '{print TOTAL,$1,$1*100/TOTAL}'
awk -vOFS='\t' '{$2=$2+$4;$3=$2+1;print $0}' ${reference}_macs_confident.txt | intersectBed -a stdin -b ${reference}_macs_control.txt | wc -l | awk -vTOTAL=$TOTAL '{print TOTAL,$1,$1*100/TOTAL}'

## Step 11: Define enriched regions
# Define regions with a confident peak in any sample as the region around the peak summit
SIZE=75 # around peak summit = 151bp ~ genomic fragment length
for pair in chip_dmel-input_dmel chip_dyak_dm3-input_dyak_dm3; do
    awk -vOFS='\t' -vSIZE=$SIZE '{s=$2+$4-SIZE;e=$2+$4+SIZE;print $1,s,e}' ${pair}_macs_confident.txt
done | sort -k1,1 -k2,2n | mergeBed -i stdin > peak_regions.txt
# For each sample and each region add the ratio of chip_read_density / input_read_density
for pair in chip_dmel-input_dmel chip_dyak_dm3-input_dyak_dm3; do
    chip=$(echo $pair | sed 's/-.*//')
    input=$(echo $pair | sed 's/.*-//')
    # Maximum chip read density for each region
    gunzip -c ${chip}.density.gz | intersectBed -a peak_regions.txt -b stdin -wao | awk '{peak=$1":"$2":"$3;if(old&&peak!=old){print max[old]+0;delete max[old]};if((!max[peak])||max[peak]<$(NF-1)){max[peak]=$(NF-1)};old=peak}END{print max[old]+0}' > tmp_${chip}
    # Maximum input read density for each region
    gunzip -c ${input}.density.gz | intersectBed -a peak_regions.txt -b stdin -wao | awk '{peak=$1":"$2":"$3;if(old&&peak!=old){print max[old]+0;delete max[old]};if((!max[peak])||max[peak]<$(NF-1)){max[peak]=$(NF-1)};old=peak}END{print max[old]+0}' > tmp_${input}
    # Ratio chip/input
    paste tmp_${chip} tmp_${input} | awk '{if($2==0){print "NA"}else{print$1/$2}}' | paste peak_regions.txt - > tmp_${pair}
    mv tmp_${pair} peak_regions.txt
    rm tmp_${chip} tmp_${input}
done

## Step 12: Data normalization
# Remove regions with no reads
awk '($4!=0&&$5!=0)' peak_regions.txt > peak_regions_no0.txt
R # Enter R
library(preprocessCore) # Load library
table_pre_norm=read.table("peak_regions_no0.txt") # Load table
table_post_norm=normalize.quantiles(as.matrix(table_pre_norm[,4:5])) # Normalize table
write.table(cbind(table_pre_norm[,1:3],signif(table_post_norm)),"peak_regions_norm.txt",quote=F,sep="\t",row.names=F,col.names=F) # Save table
q()
n

## Step 13: Quantitative changes
# Calculate log2(change)
grep -v "NA" peak_regions_norm.txt | awk -vOFS='\t' '{print $0,log($4/$5)/log(2)}' > peak_regions_norm_log2.txt
# Regions 2 fold higher in Dmel than Dyak
awk '($6>=1)' peak_regions_norm_log2.txt > peak_regions_norm_log2_decrease.txt
# Regions with no quantitative changes (within 2 fold)
awk '($6>-1&&$6<1)' peak_regions_norm_log2.txt > peak_regions_norm_log2_invariant.txt
# Regions 2 fold lower in Dmel than Dyak
awk '($6<=-1)' peak_regions_norm_log2.txt > peak_regions_norm_log2_increase.txt
# Count number of regions
wc -l peak_regions_norm_log2_*.txt
