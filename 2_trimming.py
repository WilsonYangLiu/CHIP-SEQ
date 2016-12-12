help=""" 2nd step: trim the low quality base """
import sys, string
import os
f = open(sys.argv[1], 'r')	# the list of FastQ file
allfile = f.readlines()
for i in range(len(allfile)):
	word = allfile[i].split('.')
	cmd = "gunzip -c "+allfile[i][:-1]+" | fastq_quality_filter -Q33 -v -q 20 -p 75 -o "+word[0]+"_trimed.fastq"
	os.system(cmd) 
