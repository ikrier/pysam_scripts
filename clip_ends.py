#!/usr/bin/python

import sys
import pysam
import bbcflib
from bbcflib.track import track
import clip_functions
from clip_functions import *



samfile = pysam.AlignmentFile(sys.argv[2], "rb")
pairedreads = pysam.AlignmentFile(sys.argv[3], "wb", template=samfile)


t = track(sys.argv[1],chrmeta="hg19")
t.chrmeta["chrMT"]=t.chrmeta["chrM"]
del t.chrmeta["chrM"]

for chr in t.chrmeta:
	chrom=chr[3:]
	print chrom
	my_list_starts = []
	my_list_ends= []
	for feature in t.read(str(chrom)):
		my_list_starts.append(feature[1])
		my_list_ends.append(feature[2]-1)
	#print my_list_starts, my_list_ends
	for read in samfile.fetch(str(chrom)):
		if read.is_reverse:
			if read.is_unmapped==False:
				clip_right(read) #We could maybe check if the right side is indeed an amplicon coordinate but here let's assume they all are
				if read.reference_start in my_list_starts:
					clip_left(read)
		else:
			if read.is_unmapped==False:
				clip_left(read) #idem for assuming reads start at amplicon edge
				if read.reference_end-1 in my_list_ends:
					clip_right(read)
	#	print read
		pairedreads.write(read)

pairedreads.close()
samfile.close()
