# Analysis of quantitation changes. (control sample: cell_lysate)
## Define enriched regions. Collapse all peak regions that are independently called in any of the different samples by computing the union of all peak coordinates. Score each region for each sample by the highest read counts in this region normalized to the total number of mapped reads in each sample and to number of reads at that position in the corresponding input sample(score even samples that do not have a peak in this region)
# Define regions with a confident peak in any sample as the region around the peak summit
SIZE=75 # around peak summit = 151bp ~ genomic fragment length
echo -n "Define peak regions. Time:"
date
for sample in H3K56ac H3K56HMG; do
    awk -vOFS='\t' -vSIZE=$SIZE '{s=$2+$4-SIZE;e=$2+$4+SIZE;print $1,s,e}' ${sample}_confident
done | sort -k1,1 -k2,2n | mergeBed -i stdin > peak_regions.txt
echo "Done. Next, caculate the ratio of chip and input density"
# For each sample and each region add the ratio of chip_read_density / input_read_density
for pair in H3K56ac-cell_lysate H3K56HMG-cell_lysate; do
    chip=$(echo $pair | sed 's/-.*//')
    input=$(echo $pair | sed 's/.*-//')
    # Maximum chip read density for each region
    echo -n "Maximum chip read density for each region. Time:"
    date
    gunzip -c ${chip}.density.gz | intersectBed -a peak_regions.txt -b stdin -wao | awk '{peak=$1":"$2":"$3;if(old&&peak!=old){print max[old]+0;delete max[old]};if((!max[peak])||max[peak]<$(NF-1)){max[peak]=$(NF-1)};old=peak}END{print max[old]+0}' > tmp_${chip}
    # Maximum input read density for each region
    echo -n "Maximum input read density for each region. Time:"
    date
    gunzip -c ${input}.density.gz | intersectBed -a peak_regions.txt -b stdin -wao | awk '{peak=$1":"$2":"$3;if(old&&peak!=old){print max[old]+0;delete max[old]};if((!max[peak])||max[peak]<$(NF-1)){max[peak]=$(NF-1)};old=peak}END{print max[old]+0}' > tmp_${input}
    # Ratio chip/input
    echo "Caculate the ratio"
    paste tmp_${chip} tmp_${input} | awk '{if($2==0){print "NA"}else{print $1/$2}}' | paste peak_regions.txt - > tmp_${pair}
    mv tmp_${pair} peak_regions.txt	# the last 2 columns stands for the ratio of each pair (H3K56ac-cell_lysate H3K56HMG-cell_lysate)
    #rm tmp_${chip} tmp_${input}
done

:<<BLOCK
## Data normalization. For comparisons for which a constant number of binding sites is expected in both samples, remove nonmappable region (i.e. regions without any read in one of the samples) and normalize the peak heights using quantile normalization. (Optional. Otherwise, proceed directly to next step)
# Remove regions with no reads 
awk '($4!=0&&$5!=0)' peak_regions.txt > peak_regions_no0.txt
R # Enter R
library(preprocessCore) # Load library
table_pre_norm=read.table("peak_regions_no0.txt") # Load table
table_post_norm=normalize.quantiles(as.matrix(table_pre_norm[,4:5])) # Normalize table
write.table(cbind(table_pre_norm[,1:3],signif(table_post_norm)),"peak_regions_norm.txt",quote=F,sep="\t",row.names=F,col.names=F) # Save table
q()
n

## Quantitative changes. Compute the differences between peak heights as log2 fold change.
# Calculate log2(change)
grep -v "NA" peak_regions_norm.txt | awk -vOFS='\t' '{print $0,log($4/$5)/log(2)}' > peak_regions_norm_log2.txt
# Regions 2 fold higher in H3K56ac than H3K56HMG
awk '($6>=1)' peak_regions_norm_log2.txt > peak_regions_norm_log2_decrease.txt
# Regions with no quantitative changes (within 2 fold)
awk '($6>-1&&$6<1)' peak_regions_norm_log2.txt > peak_regions_norm_log2_invariant.txt
# Regions 2 fold lower in H3K56ac than H3K56HMG
awk '($6<=-1)' peak_regions_norm_log2.txt > peak_regions_norm_log2_increase.txt
# Count number of regions
wc -l peak_regions_norm_log2_*.txt
MG-cell_lysate)
    rm tmp_${chip} tmp_${input}
done
BLOCK
