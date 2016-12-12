help = """ 3rd-1 step: align reads to the reference with bowtie, generate sam file. """
import sys, string
import os
f = open(sys.argv[1], 'r')	# the list of FastQ file
allfile = f.readlines()
for i in range(len(allfile)):
	word = allfile[i].split('.')
	cmd = "bowtie -p 6 -v 2 ./index/hg19 "+allfile[i][:-1]+" -S 2> "+word[0]+".out > "+word[0]+".sam"
	os.system(cmd)
