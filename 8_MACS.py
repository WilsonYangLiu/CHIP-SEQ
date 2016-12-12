help = """ 4th step: peek calling by using MACS """
import sys, string
import os
f = open(sys.argv[1], 'r')	# the list of sorted_rmdup_bam file, one file per line. The last line stands for the control sample
allfile = f.readlines()
for i in range(len(allfile)-1):
		word = allfile[i].split('.')
		cmd = "macs14 -t "+allfile[i][:-1]+" -c "+allfile[-1][:-1]+" --format BAM  --gsize hs --name "+word[0]+" --nomodel --shiftsize 73 --bdg --single-profile --diag &> "+word[0]+"_MACS.out"
		os.system(cmd)
