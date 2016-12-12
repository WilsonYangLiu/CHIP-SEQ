	# Genome-wide maps of H3K4me2/3 in prostate cancer cell line LNCaP
	# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20042
	# http://www2.uef.fi/documents/1698400/2466431/Macs2/f4d12870-34f9-43ef-bf0d-f5d087267602
	
    cd ~/CHIPseq_test/
    mkdir GSE20042_H3K4me2_3 && cd GSE20042_H3K4me2_3
    mkdir rawData && cd rawData
    for ((i=146;i<153;i++)) ;do wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP002/SRP002077/SRR037$i/SRR037$i.sra;done

    ls *sra |while read id; do ~/biosoft/sratoolkit/sratoolkit.2.6.3-centos_linux64/bin/fastq-dump $id;done
    rm *sra
    ls *.fastq | while read id ; do ~/biosoft/fastqc/FastQC/fastqc $id;done
    mkdir QC_results
    mv *zip *html QC_results/
    ##接下来做比对
    ## cat >run_bowtie2.sh  运行这个脚本批量做alignment
    ls *.fastq | while read id ;
    do
    echo $id
    ~/biosoft/bowtie/bowtie2-2.2.9/bowtie2 -3 5 -p 8 -x ~/biosoft/bowtie/hg19_index/hg19 -U $id   -S ${id%%.*}.sam  2>${id%%.*}.align.log;
    samtools view -bhS -q 30  ${id%%.*}.sam > ${id%%.*}.bam  ## -F 1548 https://broadinstitute.github.io/picard/explain-flags.html
    samtools sort   ${id%%.*}.bam ${id%%.*}.sorted  ## prefix for the output
    samtools index ${id%%.*}.sorted.bam
    done

	# 下载GEO的核小体定位点(peaks)结果
    tar xvf GSE20042_RAW.tar
    ls *gz |xargs gunzip
    wc -l *txt

	# 这个peakView.R代码很简单，就是用samtools depth命令提取每个peaks区域的坐标，然后画曲线即可
	Rscript ~/CHIPseq_test/peakView.R GSM503907_LNCaP_H3K4me3_DHT_4h_normalized_peak.bed  ../rawData/SRR037152.sorted.bam
	
	ls *sorted.bam |while read id;do ( nohup time ~/.local/bin/macs2 callpeak -t $id -f BAM -g hs -n ${id%%.*} 2>${id%%.*}.masc2.log &) ;done  

	## 这里批量对7个测序文件做peaks callling 
	mkdir ../MACS2results
	mv *bed *xls *Peak *r  ../MACS2results
	cd ../MACS2results
	ls *.xls | while read id ;
	do
	echo $id
	grep ‘^chr\S’ $id |perl -alne ‘{print “$F[0]\t$F[1]\t$F[2]\t$F[9]\t$F[7]\t+”}’ >${id%%.*}.bed
	done
	
	# 浏览peaks
	Rscript ~/CHIPseq_test/peakView.R SRR037152_peaks.bed  ../rawData/SRR037152.sorted.bam


###
#!/bin/bash

for i in $1/*
do
	j=`basename $i .png`
	k=`dirname $i`
	#echo "j is $j, k is $k"
	convert -contrast -thumbnail 110 "$i" $k/thumb.$j.png
done
