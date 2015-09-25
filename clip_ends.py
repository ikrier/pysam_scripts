#!/usr/bin/python
'''
The code should do the following :
	- read in a bam alignment file
	- read in a list of amplicon regions
	- clip the first 5 aligned bases of each read
	- if the end of the read is at an amplicon end, clip the last 5 aligned bases
Future improvements could include :
	- keeping deletion structures by using a placeholder N before each deleted region
	- verifying if the start of a read is at the start of an amplicon
	- veriying if reads are paired
	- allowing some uncertainty at the end of the amplicons since there's a 1/4 chance of aligning one more base, 
	  and even more if the following bases are taken into account
=======
License
=======
This code is released under the GNU General Public License 3.0. A copy
of this license is in the LICENSE.txt file.
copyright Irina Krier 2015
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import sys
import pysam
import bbcflib
from bbcflib.track import track
import clip_functions
from clip_functions import *
import argparse

parser = argparse.ArgumentParser(
    description='Clip ends of reads according to amplicon structure')
parser.add_argument('Amplicons bed file', metavar='amplicon_file.bed', type=str, nargs=1,
                   help='A bed file containing the coordinates of the amplicons used in the design')
parser.add_argument('Alignments bam', metavar='in_file.bam', type=str, nargs=1,
                   help='A bam file containing the alignments to process')
parser.add_argument('Output file', metavar='out_file.bam', type=str, nargs=1,
                   help='A bam file containing the output')

args=parser.parse_args()
myargs=vars(args)

samfile = pysam.AlignmentFile(myargs["Alignments bam"][0], "rb")
pairedreads = pysam.AlignmentFile(myargs["Output file"][0], "wb", template=samfile)

t = track(myargs["Amplicons bed file"][0],chrmeta="hg19")
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
