##########################################################################################
### 1.  Run fastqc to generate basic qc stats on sequence file                         ###
##########################################################################################
fastqc /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_Cell_Iysates/Cell_Iysates_AGTTCC_L001_R1_001.fastq --outdir=/xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_Cell_Iysates/

##########################################################################################
### 2.  Filter the sequence files for quality using fastx toolkit                      ###
##########################################################################################

fastq_quality_filter -Q33 -q 20 -p 75 -i /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_Cell_Iysates/Cell_Iysates_AGTTCC_L001_R1_001.fastq -o /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Cell_Iysates_1.fastq 
fastq_quality_filter -Q33 -q 20 -p 75 -i /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_Cell_Iysates/Cell_Iysates_AGTTCC_L002_R1_001.fastq -o /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Cell_Iysates_2.fastq 

fastq_quality_filter -Q33 -q 20 -p 75 -i /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_H3K56ac/H3K56ac_CCGTCC_L001_R1_001.fastq -o /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/H3K56ac.1.fastq
fastq_quality_filter -Q33 -q 20 -p 75 -i /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_H3K56ac/H3K56ac_CCGTCC_L002_R1_001.fastq -o /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/H3K56ac.2.fastq


fastq_quality_filter -Q33 -q 20 -p 75 -i /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_H3K56HMG/H3K56HMG_ATGTCA_L001_R1_001.fastq -o /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/H3K56HMG.1.fastq
fastq_quality_filter -Q33 -q 20 -p 75 -i /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_H3K56HMG/H3K56HMG_ATGTCA_L002_R1_001.fastq /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/H3K56HMG.2.fastq

fastq_quality_filter -Q33 -q 20 -p 75 -i /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_IgG/IgG_AGTCAA_L001_R1_001.fastq -o /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/IgG.1.fastq
fastq_quality_filter -Q33 -q 20 -p 75 -i /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/Sample_IgG/IgG_AGTCAA_L002_R1_001.fastq -o /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/IgG.2.fastq

##########################################################################################
### 3.  Map reads to the reference genome  using bowtie2                               ###
##########################################################################################
DNA_HOME=/xdfs01/genoseq/geno3/gtt/ChIP-seq/data

cat ${DNA_HOME}/sample.txt |while read line
do
bowtie2 -p 5 -x /xdfs01/genoseq/geno3/gtt/hg19/hg19 -U ${DNA_HOME}/"$line"_1.fastq  -S ${DNA_HOME}/$line.1.sam
bowtie2 -p 5 -x /xdfs01/genoseq/geno3/gtt/hg19/hg19 -U ${DNA_HOME}/"$line"_2.fastq  -S ${DNA_HOME}/$line.2.sam
done
##########################################################################################
### 4.  Use samtools to make bam, sort, rmdup and index                                ###
##########################################################################################
DNA_HOME=/xdfs01/genoseq/geno3/gtt/ChIP-seq/data

if several file,combine together:
cat replicate1.bed replicate2.bed replicate3.bed >all_replicates.bed

cat ${DNA_HOME}/sample.txt |while read line
do

wait

date
echo "!start sort and rmdup"


cd /xdfs01/genoseq/geno3/gtt/ChIP-seq/data
samtools view -bSh $line.1.sam -o $line.1.bam
samtools sort $line.1.bam $line.1.sort
samtools rmdup $line.1.sort.bam $line.1.sort.rmdup.bam
samtools index $line.1.sort.rmdup.bam

samtools view -bSh $line.2.sam -o $line.2.bam
samtools sort $line.2.bam $line.2.sort
samtools rmdup $line.2.sort.bam $line.2.sort.rmdup.bam
samtools index $line.2.sort.rmdup.bam

done

##########################################################################################
### 5.  Run MACS v1.4 to call peaks and generate coverage files                        ###
##########################################################################################

wait
cd /xdfs01/genoseq/geno3/gtt/ChIP-seq/data

macs14 -t Cell_Iysates.sort.rmdup.bam -c IgG.sort.rmdup.bam -n Cell.input -f BAM -g hs  --nomodel --shiftsize 73 -w -S

macs14 -t H3K56ac.sort.rmdup.bam -c IgG.sort.rmdup.bam -n H3K56ac.input -f BAM -g hs --nomodel --shiftsize 73 -w -S

macs14 -t H3K56HMG.sort.rmdup.bam -c IgG.sort.rmdup.bam -n H3K56HMG.input -f BAM -g hs --nomodel --shiftsize 73 -w -S
##########################################################################################
### 6.  directly upload wig files  for UCSC                                            ###   http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1%3A97072956-116072955&hgsid=422441921_pjrrvZJXFmWhYl4D42v9TCfwabh7
##########################################################################################

Upload the wig and bed files to UCSC and take a look at the data (http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1%3A118259752-120159751&hgsid=422441921_pjrrvZJXFmWhYl4D42v9TCfwabh7)

##########################################################################################
### 7.  Use bedtools to annotated peaks to hg19 features                               ###
##########################################################################################
# nearest TSS
closestBed -d -t first -a input_peaks.bed -b /xdfs01/genoseq/geno3/gtt/ChIP-seq/data/refGene_hg19_TSS.bed > input_peaks.bed.tss (cat filename1 | tr -d "\r" > newfile)
# nearest EXON
closestBed -t first -a input_peaks.bed -b refseq_exons.bed > input_peaks.bed.exon
# nearest genebody
closestBed -t first -a input_peaks.bed -b refseq_genebody.bed > input_peaks.bed.genebody
...element overlap,nearest repeat



>sed '{N;s/\n/\t/}' test.txt >test1.txt (两行变成一行)
