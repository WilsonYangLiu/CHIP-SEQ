help = """ extract the sam file headers. """
import sys, string
import os
f = open(sys.argv[1], 'r')	# the list of sam file
allfile = f.readlines()
for i in range(len(allfile)):
	word = allfile[i].split('.')
	cmd = "samtools view -H "+allfile[i][:-1]+" > "+word[0]+"_head.sam"
	os.system(cmd)
