#!/usr/bin/env python

import re
import sys

def main():

	binFile = sys.argv[1]		# say bin.5k.boundary.sorted.txt
	outputFile = sys.argv[2]	# say bin.5k.boundary.gc.txt

	dict_seq = {}
	temp_seq = []
	thisChr = 'chr1'
	with open('/home/sunjh/data3/refs/hg38/analysisSet/hg38.analysisSet.fa', 'r') as f:
		for x in f:
			if x.startswith('>'):				
				seq = ''.join(temp_seq).upper()
				dict_seq[thisChr] = seq
				temp_seq = []
				thisChr = x.strip('>').strip()
			else:
				temp_seq.append(x.rstrip())
		seq = ''.join(temp_seq).upper()
		dict_seq[thisChr] = seq
	
	#print('Part of chr1 sequence: ', dict_seq['chr1'][1:1000])
	print('Finished constructing sequence dictionary!')

	OUT = open(outputFile, 'w')
	print('chrom', 'start', 'absStart', 'end', 'length', 'mappable', 'gcContent', sep='\t', file=OUT)

	with open(binFile, 'r') as f:
		for x in f:
			x = x.strip().split('\t')
			binChr = x[0]
			binStart = int(x[1])
			binEnd = int(x[3])
			temp_seq = dict_seq[binChr][(binStart):(binEnd+1)]
			gcContent = float(len(re.findall("[CG]", temp_seq))) / float(len(re.findall("[ACGT]", temp_seq)))
			print(*x, gcContent, sep='\t', file=OUT)	

	OUT.close()


if __name__ == "__main__":
	main()
