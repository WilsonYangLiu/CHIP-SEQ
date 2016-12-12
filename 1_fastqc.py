help=""" 1st step: Quality Control using FastQC """
import sys, string
import os
f = open(sys.argv[1], 'r')	# the list of FastQ file
allfile = f.readlines()
for i in range(len(allfile)):
	word = allfile[i].split('.')
	mkdir = "mkdir "+word[0]+"_qc"
	os.system(mkdir)
	cmd = "fastqc -o ./"+word[0]+"_qc "+allfile[i][:-1]
	os.system(cmd)
