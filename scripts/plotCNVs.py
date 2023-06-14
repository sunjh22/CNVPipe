#! /usr/bin/env python

import os
from utils import readCNV2DataFrame, samplotPlot

cnvFile = snakemake.input[0]
control = snakemake.params[0]
smp_name = os.path.basename(cnvFile).split('.')[0]
cnv_df = readCNV2DataFrame(cnvFile)
for _index, row in cnv_df.iterrows():
    samplotPlot(row, [smp_name], control)

