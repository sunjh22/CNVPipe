#! /usr/bin/env python
# Find recurrent CNVs among all samples

from utils import mergeCNVFromTools
import os

def getCNVs(inputFiles):
    new_cnvs = []
    for inputFile in inputFiles:
        with open(inputFile, 'r') as f:
            sample = os.path.basename(inputFile).split('.')[0]
            for cnv in f:
                if cnv.startswith('chrom'):
                    continue
                cnv = cnv.strip().split('\t')
                chrom = cnv[0]
                start = int(cnv[1])
                end = int(cnv[2])
                cn = int(cnv[3])
                dupScore = int(cnv[6])
                toolNum = int(cnv[8])
                cnvfilter = cnv[9]
                goodScore = int(cnv[10])
                if toolNum >=3 and cnvfilter == 'True' and dupScore == 100 and goodScore == 100:
                    new_cnv = [chrom, start, end, cn, sample]
                    new_cnvs.append(new_cnv)

    return new_cnvs

def main(inputFiles, outputFile, sampleNumThe=5):
    mergedCnvs = mergeCNVFromTools(getCNVs(inputFiles))
    
    with open(outputFile, 'w') as f:
        print('chromosome', 'start', 'end', 'cn', 'samples', 'sampleNum', 'accumScore', sep='\t', file=f)
        for cnv in mergedCnvs:
            sampleNum = cnv[5]
            if sampleNum < int(sampleNumThe):
                continue
            print(*cnv, sep='\t', file=f)


if __name__ == '__main__':

    inputFiles = snakemake.input
    outputFile = snakemake.output[0]
    sampleNumThe = snakemake.params[0]  # threshold of the number of samples to keep the recurrent CNVs
    
    main(inputFiles, outputFile, sampleNumThe=sampleNumThe)