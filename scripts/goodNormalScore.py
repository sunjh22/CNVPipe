#! /usr/bin/env python

from utils import readCNVFile, calculateOverlapProp4Region

def calculateOverlapScore(target_cnv, cnv_list, reciprocal_prop=0.3):
    """
    Calculate the accumulative overlap between a target CNV and a list of cnvs, only cnvs (in cnv 
    list) with over 30% reciprocal overlap with target CNV are selected. Reduce the score 100 by the
    overlap fraction and the number of overlapped cnvs.
    - target_cnv: a single CNV with chromosome, start and end
    - cnv_list: a list of cnv regions. Could be a list of bad regions or common normal CNVs
    Return a score, with 100 means no overlap with cnv list at all.
    """
    accumLen = 0    # accumulative length of overlap region between cnv and bad
    i = 0

    for cnv in cnv_list:
        overlap, prop1, prop2 = calculateOverlapProp4Region(target_cnv, cnv)
        if min(prop1, prop2) > reciprocal_prop:
            i += 1
            accumLen += overlap

    accumProp = round(accumLen * 100 / (int(target_cnv[2]) - int(target_cnv[1])))
    score = 100 - accumProp - i * 2
    # print("This CNV totally overlaps with {:d} CNV regions\n".format(i))
    return score


inputFile = snakemake.input[0]
badListFile = snakemake.params[0]       # centromere, telomere etc. by default, we use blacklist from 10x
lowMapFile = snakemake.params[1]        # low mappable regions
normalCNVFile = snakemake.params[2]     # common SVs in normal population
outputFile = snakemake.output[0]

bad_list = readCNVFile(badListFile, tool='Bad')
lowMap_list = readCNVFile(lowMapFile, tool='Bad')
normal_list = readCNVFile(normalCNVFile, tool='Normal')

with open(inputFile, 'r') as f, open(outputFile, 'w') as g:
    print('chrom\tstart\tend\tcn\tcnv\tAS\tDS\tdhfc\tdhbfc\tdhffc\ttools\ttoolNum\tgc\tCNVfilter\tGS\tMS\tNS', file=g)
    for x in f:
        if x.startswith('chrom'):
            continue
        cnv = x.strip().split('\t')[:3]
        bad_score = calculateOverlapScore(cnv, bad_list, reciprocal_prop=0.3)
        map_score = calculateOverlapScore(cnv, lowMap_list, reciprocal_prop=0.3)
        normal_score = calculateOverlapScore(cnv, normal_list, reciprocal_prop=0.5)
        print(x.strip(), bad_score, map_score, normal_score, sep='\t', file=g)
    