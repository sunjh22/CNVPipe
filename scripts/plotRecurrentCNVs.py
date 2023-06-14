#! /usr/bin/env python

from utils import readCNV2DataFrame, samplotPlot

cnvFile = snakemake.input[0]
control = snakemake.params[0]
# Alternative: only show top10 recurrent CNVs with the most number of samples overlapped. Sort the
# DataFrame based on sample number and select top10 CNVs.
cnv_df = readCNV2DataFrame(cnvFile, recurrent=True)
for _index, row in cnv_df.iterrows():
    smp_name = row[-1].split(',')[:10]
    samplotPlot(row[:-1], smp_name, control, recurrent=True)
