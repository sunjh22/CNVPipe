#! /usr/bin/env python

# This script is used to merge consective bins in mops results

import sys

def larger2(cn1, cn2):
    if (cn1 > 2 and cn2 > 2) or (cn1 < 2 and cn2 < 2):
        return True
    else:
        return False


def merge(infile, outfile):
    with open(infile, 'r') as f, open(outfile, 'w') as g:
        _ = f.readline()
        firstLine = f.readline().strip().split('\t')
        tmpChrom, tmpStart, tmpEnd, tmpCN = firstLine[0], int(firstLine[1]), int(firstLine[2]), int(firstLine[3])
        for line in f:
            line = line.strip().split('\t')
            chrom, start, end, cn = line[0], int(line[1]), int(line[2]), int(line[3])
            if chrom == tmpChrom and larger2(tmpCN, cn) and start == tmpEnd+1:
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