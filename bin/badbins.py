#! /usr/bin/python

#set -ue

IN1 = open("hg38_bwt2_76k.bin.boundaries.sorted.txt", 'r')
IN2 = open("hg38.badRegions.txt",'r')
OUT = open("hg38_bwt2_k76.12k.badbins.v2.txt",'w')

IN2.readline()
chrname = []
badPos = []

for line in IN2.readlines():
	x = line.strip().split('\t')
	Chr = x[0]
	chrname.append(Chr)
	badPos.append(x[1:3]) # record the start and end position of bad regions

count = 0
count1 = 0
count2 = 0
for x in IN1.readlines():
	count2 += 1
	arow = x.strip().split('\t')
	binStart = arow[1]
	binEnd = arow[3]
	thischr = arow[0]
	binLen = arow[4]
	lo = True
	if thischr == 'chrX':
		thischr = 'chr23'
	elif thischr == 'chrY':
		thischr = 'chr24'
	for i in range(len(chrname)):
		if arow[0] == chrname[i]:
			if long(binStart) > long(badPos[i][0]) and long(binEnd) < long(badPos[i][1]):	# bins in bad regions
				OUT.write('1\n')
				count1 +=1
				print chrname[i]+'\t'+badPos[i][0]+'\t'+badPos[i][1]+'\t'+binStart+'\t'+binEnd
				lo = False
			elif long(binStart) < long(badPos[i][1]) and long(binEnd) > long(badPos[i][1]) and (long(badPos[i][1])-long(binStart))/long(binLen) > 0.05:	# bins overlap 70% with bad regions
				OUT.write('1\n')
				count1 +=1
				print chrname[i]+'\t'+badPos[i][0]+'\t'+badPos[i][1]+'\t'+binStart+'\t'+binEnd
				lo = False
			elif long(binStart) < long(badPos[i][0]) and long(binEnd) > long(badPos[i][0]) and (long(binEnd)-long(badPos[i][0]))/long(binLen) > 0.05:	#bins overlap 70% with bad regions
				OUT.write('1\n')
				count1 +=1
				print chrname[i]+'\t'+badPos[i][0]+'\t'+badPos[i][1]+'\t'+binStart+'\t'+binEnd
				lo = False
			else:
				continue
	if lo is True:
		count += 1
		OUT.write('0\n')
	else:
		continue


print count2
print count
print count1
IN1.close()
IN2.close()
OUT.close()
