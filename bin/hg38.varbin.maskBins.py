#!/usr/bin/env python

import sys

if len(sys.argv) != 4:
	sys.exit("Argument: hg38.varbin.bins.gc.content.file gaps bad.bins.output")

bins = sys.argv[1]	# say hg38.varbin.5k.gc.content.samtools.k150.sort.txt
gaps = sys.argv[2]	# say GRCh38.centromere.telomere.extend.c_5kb.t_1kb.txt
output = sys.argv[3]    # say hg38.varbin.5k.bad.bins.txt OR GRCh38.badRegions.txt
cri = 0.10

chr2gaps = {}
GAP = open(gaps, 'r')

for line in GAP:
	if line.startswith('#'):
		continue
	arow = line.strip().split('\t')
#	Chr = line.split('\t')[0].upper().replace('CHR','')
	Chr = arow[0].upper().replace('CHR','')
	Chr = Chr.upper().replace('CHRX','23')
	Chr = Chr.upper().replace('CHRY','24')
	start = int(arow[1])
	end = int(arow[2])
	try:
		chr2gaps[Chr].append([start,end])
	except KeyError:
		chr2gaps[Chr] = [[start,end]]
GAP.close()

mask = []
BIN = open(bins, 'r')

for line in BIN:
	if line.startswith('bin'):
		continue
	arow = line.strip().split('\t')
	Chr = arow[0].upper().replace('CHR','')
	Chr = Chr.upper().replace('CHRX','23')
	Chr = Chr.upper().replace('CHRY','24')
	start = int(arow[1]) + 1
	end = int(arow[3]) + 1
	intersect = False
	try:
		for gap in chr2gaps[Chr]:
			max_start = max(gap[0], start)
			min_end = min(gap[1], end)
			overlap = min_end - max_start + 1
			if float(overlap) / (end - start + 1) > cri:
				intersect = True
				break
	except KeyError:
		pass
	if intersect:
		mask.append('1\n')
	else:
		mask.append('0\n')
BIN.close()
MASK = open(output, 'w')
MASK.writelines(mask)
MASK.close()

