#!/usr/bin/env python3
"""
Written by Will Overholt
2017-10-09

Very crudely hacked from the QIIME 1.9.1 script library
http://qiime.org/
"""

input_map = snakemake.input.full_map
output_biom = snakemake.output.biom
output_table = snakemake.output.table


import sys
from collections import defaultdict
from datetime import datetime
from biom.table import Table

#Parse OTU Map
def parse_otu_map(otu_map_fp, delim="_"):
    result = defaultdict(int)
    sample_ids = []
    sample_id_idx = {}
    otu_ids = []
    otu_count = 0
    sample_count = 0
    with open(otu_map_fp, "r") as otu_map_fp:
        for line in otu_map_fp:
            fields = line.strip().split('\t')
            otu_id = fields[0]
            for seq_id in fields[1:]:
                sample_id = seq_id.split("_")[0]
                try:
                    sample_index = sample_id_idx[sample_id]
                except KeyError:
                    sample_index = sample_count
                    sample_id_idx[sample_id] = sample_index
                    sample_count += 1
                    sample_ids.append(sample_id)
                # {(row,col):val}
                result[(otu_count, sample_index)] += 1
            otu_count += 1
            otu_ids.append(otu_id)
        return result, sample_ids, otu_ids

#Generate the table structure for the biom package
def make_otu_table(otu_map_fp, table_id=None):
    data, sample_ids, otu_ids = parse_otu_map(otu_map_fp, delim="_")
    try:
        return Table(data, otu_ids, sample_ids, table_id=table_id,
                     generated_by="Will_Overholt",
                     create_date=datetime.now().isoformat())

    except ValueError as e:
        raise ValueError("Couldn't create OTU table. Is your OTU map empty?"
                         " Original error message: %s" % (str(e)))

biom_table = make_otu_table(input_map)

#Writeout the 2 file formats
with open(output_biom, 'w') as biom_file:
    #hardcoded for json instead of hd5 (can change)
    biom_table.type = "OTU table"

    biom_table.to_json("Will_Overholt", biom_file)
#This is just the old format OTU table, can always convert it using the biom commands
print(biom_table, file = open(output_table, "w"))
