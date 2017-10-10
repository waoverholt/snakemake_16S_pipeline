#!/usr/bin/env python3
"""
Written by Will Overholt
6-25-2016
"""

import sys,re,os
from collections import OrderedDict

UCFileName = snakemake.input[0]
OUTFILE = snakemake.output[0]
OTU_DICT = OrderedDict()
with open(UCFileName) as UC:
	for line in UC:
		line = line.rstrip()
		Field = line.split("\t")
		if Field[0] == 'S':
			#this is a seed sequence for the dict key
			OTU_DICT.setdefault(Field[8])
		elif Field[0] == 'H':
			if Field[9] in OTU_DICT:
				try:
					OTU_DICT[Field[9]].append(Field[8])
				except AttributeError:
					OTU_DICT[Field[9]] = [Field[8]]
			#this is a hit against a seed
			#OTU_DICT.setdefault(Field[9]).append(Field[8])
			#OTU_DICT[Field[9]].append(Field[8])

with open(OUTFILE, "w") as OUT:
	for i,OTU in enumerate(OTU_DICT, 1):
		LABEL = snakemake.params.derep_otu_label
		#i is the OTU number (number of seeds)
		#OTU is the dict key name (the name of the seed sequence)
		OUT.write("{0}{1}\t{2}".format(LABEL, i, OTU))
		if OTU_DICT[OTU]:
			#if the seed has other identical OTUs
			for s in OTU_DICT[OTU]:
				#loop through the other OTU & print them tab-sep on the same line
				OUT.write("\t" + s)
			#After looping through all values, add new line
			OUT.write("\n")
		else:
			OUT.write("\n")
