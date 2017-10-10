#!/usr/bin/env python3

import sys,re,os
from collections import OrderedDict

first_map = snakemake.input.derep_map
second_map = snakemake.input.otu_map
out_map = snakemake.output[0]

first_dict = OrderedDict()
sec_dict = OrderedDict()

with open(first_map, "r") as FM:
    for line in FM:
        line = line.rstrip()
        cluster_name = line.split("\t")
        first_dict.setdefault(cluster_name[0],[])
        seqs = cluster_name[1].split(" ")
        first_dict[cluster_name[0]] = [seq for seq in seqs]

with open(second_map, "r") as SM:
    for line in SM:
        line = line.rstrip()
        OTU_name = line.split("\t")
        sec_dict.setdefault(OTU_name[0],[])
        seqs = OTU_name[1].split(" ")
        sec_dict[OTU_name[0]] = [seq for seq in seqs]


with open(out_map, "w") as out_fp:
    for OTU in sec_dict:
        out_fp.write("{}\t".format(OTU))
        for cluster in sec_dict[OTU]:
            if cluster in first_dict:
                #print(cluster)
                for seq in first_dict[cluster]:
                    out_fp.write("{}\t".format(seq))
        out_fp.write("\n")

"""
for OTU in sec_dict:
    sys.stdout.write("{}\t".format(OTU))
    for cluster in sec_dict[OTU]:
        if cluster in first_dict:
            for seq in first_dict[cluster]:
                sys.stdout.write("{}\t".format(seq))
    sys.stdout.write("\n")
"""
