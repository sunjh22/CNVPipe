#! /usr/bin/env python

# Construct bin boundary file, keep each bin has equal mappable sites

import sys
import bisect

def main():

	bincount = 600000		# 600k bin count means resolution is around 5k
	MAP = sys.argv[1]		# say k100.mappable.bed
	GOOD = sys.argv[2]		# say k100.umap.bed
	CHROMLEN = sys.argv[3]	# say hg38.chrom.sizes.txt
	outputfile = sys.argv[4]	# say bin.5k.boundaries.txt

	# get absolute position of chromosome end
	chromlen = dict()
	with open(CHROMLEN, 'r') as f:
		for x in f:
			arow = x.rstrip().split("\t")	
			thisChrom = arow[0]
			thisChromlen = int(arow[1])
			thisAbspos = int(arow[2])
			chromlen[thisChrom] = [thisChromlen, thisAbspos]

	# get size of mappable regions in each chromosome
	chroms = dict()
	totalLength = int(0)
	with open(MAP, 'r') as f:
		for x in f:
			if x.startswith('chrM'):
				continue
			arow = x.rstrip().split("\t")
			thisChrom = arow[0]
			thisLength = int(arow[1])
			totalLength += thisLength
			chroms[thisChrom] = thisLength

	# assign different number of bins to each chromosome according to specific bincount
	chromarray = []
	bincountUsed = 0
	for k, v in chroms.items():
		chromBincount = float(bincount) * (float(v) / float(totalLength))
		# bin number in each chromosome
		i = int(chromBincount)
		bincountUsed += i
		# r is a decimal
		r = chromBincount - i
		# k,i,r is chrom name, bincount, decimal, i.e. [['chr1', 47034, 0.72],['chr2', 50041, 0.88]]
		chromarray.append([k, i, r])

	# sort chromosomes based on decimal, the larger the decimal, the higher chance for that chromosome to get an extra bin
	a = []
	for i in chromarray:
		bisect.insort(a, (-i[2], i))	# a: [(-0.1, ['chr1', 5, 0.1]), (0.2, ['chr2', 4, -0.2])]

	chromarray = []
	for j in a:
		chromarray.append(j[1]) ##  chromarray: [['chr1', 5, 0.1], ['chr2', 4, -0.2]]. These two steps are Totally unreasonable.

	remain = bincount - bincountUsed
	#print('Remain bins', remain)
	for i in range(remain):
		chromarray[i][1] += 1

    # calculate the bin number and their length in each chromosome
	chroms2 = dict()
	for i in range(len(chromarray)):
		thischr = chromarray[i][0]
		chromlength = chroms[thischr]
		chrombins = chromarray[i][1]
		binlength = float(chromlength) / float(chrombins)
		chroms2[thischr] = [chrombins, binlength]

	#print('Bin count and length in chromosome: ', chroms2)

	OUT = open(outputfile, 'w')
	with open(GOOD, 'r') as f:
		flagChr = 'chr1'
		binCount = 0		
		currentStart = 0
		currentLength = 0
		binStart = 0
		currentExcess = 0
		for x in f:
			if x.startswith('track'):
				continue
			x = x.rstrip().split('\t')

			# each umap bed region
			thisChr = x[0]
			thisStart = int(x[1])
			thisEnd = int(x[2])
			thisSize = thisEnd - thisStart
			currentLength += thisSize

			# required bin number in current chromosome and specified bin length
			chromBins = int(chroms2[flagChr][0])
			binLength = int(chroms2[flagChr][1])
			excess = chroms2[thisChr][1] - binLength
			
			if thisChr != flagChr:
				print('Total bins in ', flagChr, ' is: ', chromBins)
				flagChr = thisChr
				binCount = 0				
				currentStart = 0
				currentLength = thisSize
				binStart = 0
				currentExcess = 0
				
			while currentLength > binLength:
				# every two or three bins, we need to shift binLength by 1 to compensate what we have ignored for decimal of binlength.
				# if we do not shift, when the resolution is high, it will greatly affect the final bin division.
				temp_binLength = binLength
				if currentExcess > 1:
					temp_binLength = binLength + 1
					currentExcess -= 1
				currentLength -= temp_binLength
				currentStart = thisEnd - currentLength
				binEnd = currentStart
				binSize = binEnd - binStart
				binStartAbspos = binStart + chromlen[thisChr][1]
				binCount += 1
				currentExcess += excess
				if binCount == chromBins:
					print('We have generated ', binCount, ' bins in this chromosome.')
					print(thisChr, binStart, binStartAbspos, chromlen[flagChr][0], binSize, temp_binLength, sep='\t', file=OUT)
					break
				else:
					print(thisChr, binStart, binStartAbspos, binEnd, binSize, temp_binLength, sep='\t', file=OUT)
				binStart = currentStart
				

	print('Bins in chrY is: ', chroms2['chrY'][0])

	OUT.close()
				

if __name__ == "__main__":
	main()
