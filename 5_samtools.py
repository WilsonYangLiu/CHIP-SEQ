help = """ 3rd-2 step: do convert,sort and index for these sam. """
import sys, string
import os
f = open(sys.argv[1], 'r')	# the list of sam file
allfile = f.readlines()
for i in range(len(allfile)):
	word = allfile[i].split('.')
	cmd1 = 'samtools view -bS '+allfile[i][:-1]+' > '+word[0]+'.bam'
	os.system(cmd1)
	cmd2 = 'samtools sort '+word[0]+'.bam '+word[0]+'_sorted'
	os.system(cmd2)
	#cmd3 = 'samtools rmdup -s '+word[0]+'_sorted.bam '+word[0]+'_sorted_rmdup.bam'
	#os.system(cmd3)
	#cmd4 = 'samtools index '+word[0]+'_sorted_rmdup.bam '+word[0]+'_sorted_rmdup.bai'
	#os.system(cmd4)
