# Restrict the dataset to the summit of the peaks +/- 100bp
perl -lane '$start=$F[1]-100 ; $end = $F[2]+100 ; print "$F[0]\t$start\t$end"' macs14_summits.bed > macs14_summits+-100.bed
# Extract the sequences for this BED file
bedtools getfasta -fi Escherichia_coli_K_12_MG1655.fasta -bed macs14_summits+-100.bed -fo macs14_summits+-100.fa
