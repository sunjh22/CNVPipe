#! /usr/bin/env python
# This script is wrote to integrate single cell CNA profiles
# in Snakemake pipeline

from collections import defaultdict

inputfile = snakemake.input[:]
outputfile = snakemake.output[0]
medianFile = snakemake.params[0]

samples = []
with open(medianFile, 'r') as f:
	for x in f:
		x = x.strip().split('\t')
		sample = x[0]
		totalReads = float(x[1])
		medianCount = float(x[2])
		if totalReads > 700000 and medianCount > 35:
			samples.append(sample)

samples = set(samples)

chrInfo = defaultdict(list)
chrInfo['chr'].append('chr')
chrInfo['pos'].append('pos')

matrix = []
flag = 0

for x in inputfile:
	node = x.split('/')[-1].split('.')[0]
	if node not in samples:
		print("Sample {:s} does not have good quality".format(node))
		continue
	print("Sample processing:", node)
	temp = [node]
	with open(x, 'r') as f:
		for line in f:
			# remove CNAs in chromosome Y
			if line.startswith('chrom') or line.startswith('24'):
				continue
			line = line.strip().split('\t')
			if flag == 0:
				chrInfo['chr'].append(line[0])
				chrInfo['pos'].append(line[1])
			integerCN = int(line[8])
			temp.append(integerCN)
	if chrInfo['chr'] not in matrix:
		matrix.append(chrInfo['chr'])
		matrix.append(chrInfo['pos'])
	matrix.append(temp)
	flag = 1

t_matrix = list(zip(*matrix))

with open(outputfile, 'w') as f:
	for i in range(len(t_matrix)):
		print(*t_matrix[i], sep = '\t', file=f)

print('Totally {:d} samples are integrated'.format(len(samples)))