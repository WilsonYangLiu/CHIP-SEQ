# !bin/bash

# FastQC
## change to the working directory
path=/export/home/liuwx/Desktop/Project_1832724646/
cd $path

# change into every child directory and run pipline
:<<BLOCK
   i'am annotation!
BLOCK
filename=(`ls`)
for var in ${filename[@]};do
wkdir="$path""$var"
echo $wkdir
cd $wkdir
## use fastqc to get QC scores
#fastqc *.fastq.gz
## use fastq_quality_filter tool to trim the read
fastq_quality_filter -q 22 -p 75 -i 
done





