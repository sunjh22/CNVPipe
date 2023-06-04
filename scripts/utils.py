#! /usr/bin/env python

# Includes all functions used in CNVPipe

def testSameCNVType(cn1, cn2):
    """Test if two CNVs are the same type"""
    if (cn1 > 2 and cn2 > 2) or (cn1 < 2 and cn2 < 2):
        return True
    else:
        return False


def readCNVFile(infile, tool):
    """Read CNVKit output file into a list"""
    cnvList = []
    with open(infile) as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            if tool == 'CNVKit':
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[5])
            elif tool in ['MOPS', 'merge']:
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
            elif tool in ['Smoove', 'Delly']:
                if line.find('_') != -1 and not line.startswith('NC'):
                    continue
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
            elif tool in ['Bad', 'Normal']:
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), 2
            else:
                raise ValueError('No tool is indicated, please indicate one!')
            if cn == 2:
                continue
            cnvList.append([chrom, start, end, cn])

    return cnvList


def testOverlapped(cnv1, cnv2):
    """ Test whether two CNVs are overlapped
    >>> testOverlapped(['chr1',1,10], ['chr1',10,20])
    False

    >>> testOverlapped(['chr1',5,15], ['chr1',10,20])
    True

    >>> testOverlapped(['chr1',12,20], ['chr1',10,20])
    True

    >>> testOverlapped(['chr1',15,25], ['chr1',10,20])
    True

    >>> testOverlapped(['chr1',20,30], ['chr1',10,20])
    False

    >>> testOverlapped(['chr1',8,25], ['chr1',10,20])
    True

    >>> testOverlapped(['chr2',10,20], ['chr1',10,20])
    False
    """
    c1, s1, e1 = cnv1[0], int(cnv1[1]), int(cnv1[2])
    c2, s2, e2 = cnv2[0], int(cnv2[1]), int(cnv2[2])
    cnvLen1, cnvLen2 = e1 - s1, e2 - s2
    assert cnvLen1 > 0 and cnvLen2 > 0, "The length of CNV is less than 0, please check your input CNVs."
    if c1 == c2:
        if s2 < s1 < e2 or s2 < e1 < e2 or s1 < s2 < e2 < e1:
            return True
        
    return False


def calculateOverlapProp(cnv1, cnv2):
    """ Calculate the overlapped proportion between two CNVs
    >>> calculateOverlapProp(['chr1',100,1000,1], ['chr1',500,1000,0])
    (500, 0.56, 1.0)
    """
    c1, s1, e1, cn1 = cnv1[0], int(cnv1[1]), int(cnv1[2]), int(cnv1[3])
    c2, s2, e2, cn2 = cnv2[0], int(cnv2[1]), int(cnv2[2]), int(cnv2[3])
    cnvLen1, cnvLen2 = e1 - s1, e2 - s2
    overlap = 0
    assert cnvLen1 > 0 and cnvLen2 > 0, "The length of CNV is less than 0, please check input CNVs"
    if c1 == c2:
        if testSameCNVType(cn1, cn2):
            if s2 <= s1 < e1 <= e2:
                overlap = e1 - s1
            elif s2 <= s1 <= e2 < e1:
                overlap = e2 - s1
            elif s1 <= s2 < e2 <= e1:
                overlap = e2 - s2
            elif s1 < s2 <= e1 <= e2:
                overlap = e1 - s2
            else:
                pass
    
    prop1 = round((overlap/cnvLen1), 2)
    prop2 = round((overlap/cnvLen2), 2)
    return overlap, prop1, prop2


def calculateOverlapProp4Region(cnv1, cnv2):
    """ Calculate the overlapped proportion between CNV and genomic region
    >>> calculateOverlapProp(['chr1',100,1000], ['chr1',500,1000])
    (500, 0.56, 1.0)
    """
    c1, s1, e1 = cnv1[0], int(cnv1[1]), int(cnv1[2])
    c2, s2, e2 = cnv2[0], int(cnv2[1]), int(cnv2[2])
    cnvLen1, cnvLen2 = e1 - s1, e2 - s2
    overlap = 0
    assert cnvLen1 > 0 and cnvLen2 > 0, "The length of CNV is less than 0, please check input CNVs"
    if c1 == c2:
        if s2 <= s1 < e1 <= e2:
            overlap = e1 - s1
        elif s2 <= s1 <= e2 < e1:
            overlap = e2 - s1
        elif s1 <= s2 < e2 <= e1:
            overlap = e2 - s2
        elif s1 < s2 <= e1 <= e2:
            overlap = e1 - s2
        else:
            pass
    
    prop1 = round((overlap/cnvLen1), 2)
    prop2 = round((overlap/cnvLen2), 2)
    return overlap, prop1, prop2


def mergeConsecutiveSegments(cnvList, shift=0):
    """ Merge consecutive segments and average copy number
    >>> mergeConsecutiveSegments([['chr1',1,10,1], ['chr1',10,20,0], ['chr2',1,5,3], ['chr2',5,10,4], ['chr2',10,20,1], ['chr3',1,10,4], ['chr3',10,20,8]], shift=0)
    [['chr1', 1, 20, 0], ['chr2', 1, 10, 4], ['chr2', 10, 20, 1], ['chr3', 1, 20, 6]]

    >>> mergeConsecutiveSegments([['chr1',1,10,1], ['chr1',11,20,0], ['chr2',1,5,3], ['chr2',6,10,4], ['chr2',11,20,1], ['chr3',1,10,4], ['chr3',11,20,8]], shift=1)
    [['chr1', 1, 20, 0], ['chr2', 1, 10, 4], ['chr2', 11, 20, 1], ['chr3', 1, 20, 6]]
    """
    new_cvnList = []
    tmpChrom, tmpStart, tmpEnd, tmpCN = cnvList[0]
    for i in range(1, len(cnvList)):
        chrom, start, end, cn = cnvList[i]
        if chrom == tmpChrom and testSameCNVType(tmpCN, cn) and start == tmpEnd + shift:
            tmpEnd = end
            tmpCN = round((tmpCN+cn)/2)
        else:
            new_cvnList.append([tmpChrom, tmpStart, tmpEnd, tmpCN])
            tmpChrom, tmpStart, tmpEnd, tmpCN = chrom, start, end, cn
    
    new_cvnList.append([tmpChrom, tmpStart, tmpEnd, tmpCN])
    return new_cvnList


def resolveConflictCNVs(cnvList):
    """ Find and remove conflicted CNVs
    >>> resolveConflictCNVs([['chr1',100,1000,1], ['chr1',500,1200,3], ['chr1',800,2000,4], ['chr2',100,1000,3], ['chr3',100,1000,1], ['chr3',1000,2000,4], ['chr3',1500,2500,4]])
    [['chr2', 100, 1000, 3], ['chr3', 100, 1000, 1]]
    """
    new_cnvList = []
    for x in cnvList:
        cnvLen = int(x[2]) - int(x[1])
        flag = 0
        for y in cnvList:
            if x != y:
                if testOverlapped(x[:3], y[:3]) > 0:
                    flag = 1

        # the length of identified CNV should be larger than 50bp
        if flag == 0 and cnvLen >= 100:
            new_cnvList.append(x)

    return new_cnvList


def mergeCNVFromTools(cnvList, min_threshold=0.75, max_threshold=0.95):
    """ Merge CNV result from different tools
    >>> mergeCNVFromTools([['chr1',0,1000,1,'smoove'], ['chr1',2000,4000,3,'smoove'], ['chr1',5000,6000,4,'smoove'], \
                           ['chr1',100,1000,0,'delly'], ['chr1',5500,6500,4,'delly'], ['chr2',0,1000,1,'delly'], \
                           ['chr1',200,1100,1,'cnvkit'], ['chr1',5100,6000,4,'cnvkit'], ['chr2',100,1100,1,'cnvkit'], \
                           ['chr1',200,1200,1,'cnvpytor'], ['chr1',5500,6400,4,'cnvpytor'], ['chr3',0,1100,1,'cnvpytor'], \
                           ['chr1',1000,1100,1,'mops'], ['chr1',5300,6300,4,'mops'], ['chr4',0,1000,1,'mops']])
    [['chr1', 0, 1200, 1, 'smoove,delly,cnvkit,cnvpytor,mops', 5, 225.0], ['chr1', 2000, 4000, 3, 'smoove', 1, 0.0], ['chr1', 5000, 6000, 4, 'smoove,cnvkit', 2, 90.0], ['chr1', 5300, 6500, 4, 'delly,cnvpytor,mops', 3, 141.7], ['chr2', 0, 1100, 1, 'delly,cnvkit', 2, 81.8], ['chr3', 0, 1100, 1, 'cnvpytor', 1, 0.0], ['chr4', 0, 1000, 1, 'mops', 1, 0.0]]
    """
    cnvs2 = cnvList[:]
    mergedCnvs = []
    while cnvList:
        cnv1 = cnvList.pop(0)
        cnvs2.pop(0)
        cnvs3 = cnvs2[:]
        count = 0
        accumLen = 0
        tmpCN = [cnv1[3]]
        for i, cnv2 in enumerate(cnvs3):
            # if overlap proportion is larger than the threshold, extend breakpoints
            overlap, prop1, prop2 = calculateOverlapProp(cnv1[:4], cnv2[:4])
            if min(prop1, prop2) > min_threshold or max(prop1, prop2) > max_threshold:
                # assert cnv1[-1] != cnv2[-1], f"Overlapped CNV from same tool {cnv1[-1]:s}! Please make sure these conflicts are solved before merging"
                cnv1[1:3] = [min(cnv1[1], cnv2[1]), max(cnv1[2], cnv2[2])]
                cnv1[-1] = ','.join([cnv1[-1], cnv2[-1]])
                cnvs2.pop(i-count)
                accumLen += overlap
                count += 1
                tmpCN.append(cnv2[3])

        tools = set(cnv1[-1].split(','))
        cnv1[-1] = ",".join(tools)
        cnv1.append(len(tools))
        cnv1[3] = (max(tmpCN, key=tmpCN.count))
        accumFold = round(accumLen * 100 / (int(cnv1[2]) - int(cnv1[1])), 1)   # percentage
        cnv1.append(accumFold)
        mergedCnvs.append(cnv1)
        cnvList = cnvs2[:]

    return mergedCnvs


if __name__ == '__main__':

    import doctest
    doctest.testmod()