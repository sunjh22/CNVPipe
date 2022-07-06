#! /usr/bin/env python
# Sum up the length of mappable regions in every chromosome

import sys

def main():
	
	inputfile = sys.argv[1]		# say k100.umap.bed
	outputfile = sys.argv[2]	# say k100.mappable.bed

	prevChrom = "chr1"
	prevChromMappable = 0

	with open(inputfile, 'r') as f:
		with open(outputfile, 'w') as g:
			for x in f:
				if x.startswith('track'):
					continue
				arow = x.rstrip().split("\t")
				thisChrom = arow[0]
				thisLength = int(arow[2]) - int(arow[1])

				if thisChrom == prevChrom:
					prevChromMappable += thisLength
				else:
					print(prevChrom, prevChromMappable, sep='\t', file=g)
					prevChromMappable = thisLength

				prevChrom = thisChrom

			print(prevChrom, prevChromMappable, sep='\t', file=g)


if __name__ == "__main__":
	main()
