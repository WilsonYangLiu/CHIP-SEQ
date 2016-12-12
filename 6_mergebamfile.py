help = """ merge a list of bam files, then do sort, rmdup, index for the merge file. """
import sys, string
import os
header = sys.argv[1]		# the header file
lists = sys.argv[2]		# the name of the file which cantains list of input BAM files, one file per line
######################################################################
## output the merge file with the name of the lists                ###
######################################################################
cmd1 = "samtools merge -h "+header+" -b "+lists+" "+lists+".bam"
os.system(cmd1)
######################################################################
## do sort, rmdup, index for the merge file                        ###
######################################################################
cmd2 = "samtools sort "+lists+".bam "+lists+"_sorted"
os.system(cmd2)
cmd3 = "samtools rmdup -s "+lists+"_sorted.bam "+lists+"_sorted_rmdup.bam"
os.system(cmd3)
cmd4 = "samtools index "+lists+"_sorted_rmdup.bam "+lists+"_sorted_rmdup.bai"
os.system(cmd4)
