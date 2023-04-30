#! /usr/bin/env python

# This script is used to merge consective bins in mops results

import sys

def larger2(cn1, cn2):
    if (cn1 > 2 and cn2 > 2) or (cn1 < 2 and cn2 < 2):
        return True
    else:
        return False


def readFile(infile):
    cnvList = []
    with open(infile) as f:
        for x in f:
            if x.startswith('chrom'):
                continue
            x = x.strip().split('\t')
            chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[5])
            if cn == 2:
                continue
            cnvList.append([chrom, start, end, cn])

    return cnvList

def merge(infile, outfile):
    with open(outfile, 'w') as g:
        cnvList = readFile(infile)
        tmpChrom, tmpStart, tmpEnd, tmpCN = cnvList[0]
        for i in range(1, len(cnvList)):
            cnv = cnvList[i]
            chrom, start, end, cn = cnv
            if chrom == tmpChrom and larger2(tmpCN, cn) and start == tmpEnd:
                tmpEnd = end
                tmpCN = round((tmpCN+cn)/2)
            else:
                print(tmpChrom, tmpStart, tmpEnd, tmpCN, sep='\t', file=g)
                tmpChrom, tmpStart, tmpEnd, tmpCN = chrom, start, end, cn
        
        print(tmpChrom, tmpStart, tmpEnd, tmpCN, sep='\t', file=g)
        

if __name__ == '__main__':

    inputFile = sys.argv[1]
    outputFile = sys.argv[2]

    merge(inputFile, outputFile)