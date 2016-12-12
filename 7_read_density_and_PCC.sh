# Read density visualization
# Visualize the read density with BigWig files (compressed binary version of WIG files) by extending the reads to the average length of the genomic fragments known a priori or determined during peak calling and counting the number of reads at each position in the genome normalized to the total number of mapped reads in the library.
:<<BLOCK
date
for sample in H3K56ac H3K56HMG IgG cell_lysate; do
    EXTEND=150
    GENOME=hg19.genome
    # Number of reads
    librarySize=$(samtools idxstats ${sample}.bam | awk '{total+=$3} END{print total}')
    echo "Library size: "$librarySize
    # Create density file: extend reads, calculate read density at each position and normalize the library size to 1 million reads
    echo -n "Create density file; "
    date
    bamToBed -i ${sample}.bam | awk -vCHROM=$GENOME -vEXTEND=$EXTEND -vOFS='\t' 'BEGIN{while(getline<CHROM){chromSize[$1]=$2}}{chrom=$1;start=$2;end=$3;strand=$6;if(strand=="+"){end=start+EXTEND;if(end>chromSize[chrom]){end=chromSize[chrom]}};if(strand=="-"){start=end-EXTEND;if(start<1){start=1}};print chrom,start,end}' | sort -k1,1 -k2,2n | genomeCoverageBed -i stdin -g $GENOME -d | awk -vOFS='\t' -vSIZE=$librarySize '{print $1,$2,$2+1,$3*1000000/SIZE}' | gzip > ${sample}.density.gz
    # Create WIG file
    echo -n "Create WIG file; "
    date
    gunzip -c ${sample}.density.gz | awk -vOFS='\t' '($4!=0){if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | gzip > ${sample}.wig.gz
    # Create BigWig file
    echo -n "Create BigWig file; "
    date
    wigToBigWig ${sample}.wig.gz $GENOME ${sample}.bw
    # Remove intermediate file
    rm ${sample}.wig.gz
    date
done
BLOCK

# Assessing global reproducibility and similarity
## PCC(Pearson correlation coefficient). Calculate the PCC between the normalized extended read counts at each position in the reference genome for every pair of sample
### Paired sample: H3K56ac-cell_lysate H3K56HMG-cell_lysate H3K56ac-H3K56HMG H3K56ac-IgG H3K56HMG-IgG; 
echo -e "pair of samples\tcorrelation"
for pair in H3K56ac-cell_lysate H3K56HMG-cell_lysate H3K56ac-H3K56HMG H3K56ac-IgG H3K56HMG-IgG; do
    echo -en $pair"\t"
    chip=$(echo $pair | sed 's/-.*//')
    input=$(echo $pair | sed 's/.*-//')
    paste <(gunzip -c ${chip}.density.gz) <(gunzip -c ${input}.density.gz) | awk '{if($2!=$6){exit 1};if($4!=0||$8!=0){print $4,$8}}' | ./correlation.awk
done
