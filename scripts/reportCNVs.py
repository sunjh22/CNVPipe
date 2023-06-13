#! /usr/bin/env python

import pandas as pd
import numpy as np
import sys
import os
import subprocess
from PIL import Image


def samplot_plots(cnvFile, bamFile, ctrsmps, outFile):
    """Use samplot to plot reads distribution for top10 cnvs with the highest overalll score"""
    smp_name = os.path.basename(cnvFile).split('.')[0]
    cnv_df = read_cnv(cnvFile)
    for _index, row in cnv_df.iterrows():
        samplot_plot(row, smp_name, bamFile, ctrsmps)
    
    concat_pngs(outFile)


def read_cnv(infile):
    """Read output from CNVPipe, convert into a DataFrame"""
    cnv_df = pd.read_csv(infile, sep='\t', usecols=[0,1,2,3],
                         dtype={'chrom':str, 'start':int, 'end':int, 'cn':int})
    cnv_df['cnv'] = ['DUP' if cn_value > 2 else 'DEL' for cn_value in cnv_df['cn']]
    return cnv_df.iloc[:3].drop('cn', axis=1)


def samplot_plot(cnvs, smp_name, bamFile, ctrsmps):
    chrom, start, end, cnv = cnvs
    tgt_dir = 'res/report/'+smp_name
    if not os.path.exists(tgt_dir):
        os.makedirs(tgt_dir, exist_ok=True)
    outFile = os.path.join(tgt_dir, chrom+'_'+str(start)+'_'+str(end)+'_'+cnv+'.png')
    ctr_smp_name = ' '.join(ctrsmps)
    ctrbams = ['mapped/'+sample+'.bam' for sample in ctrsmps]
    ctr_bam_name = ' '.join(ctrbams)
    # print(smp_name, ctr_smp_name, bamFile, ctr_bam_name)
    subprocess.run(['samplot', 'plot',
                    '-n', smp_name, ctr_smp_name,
                    '-b', bamFile, ctr_bam_name,
                    '-c', chrom,
                    '-s', str(start),
                    '-e', str(end),
                    '-t', cnv,
                    '-w 300',
                    '--same_yaxis_scales', 
                    '--max_coverage_points 1000',
                    '--coverage_tracktype superimpose', 
                    '--coverage-only', 
                    '-o', outFile])


def concat_pngs(outFile):
    smp_name = os.path.basename(outFile).split('.')[0]
    files = os.listdir('res/report/'+smp_name)
    img_open_list = []
    for file in files:
        img_open = Image.open('res/report/'+smp_name+'/'+file)
        if img_open.mode != 'RGB':
            img_open = img_open.convert('RGB')
        img_open_list.append(img_open)
    img1 = img_open_list[0]
    img_open_list = img_open_list[1:]
    img1.save(outFile, 'PDF', resolution=300, save_all=True, append_images=img_open_list)


if __name__ == '__main__':

    cnvFile = sys.argv[1]
    snpFile = sys.argv[2]
    bamFile = sys.argv[3]
    outFile = sys.argv[4]
    controlSamples = sys.argv[5:]

    samplot_plots(cnvFile, bamFile, controlSamples, outFile)
