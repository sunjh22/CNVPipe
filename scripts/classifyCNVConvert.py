#! /usr/bin/env python

from collections import defaultdict

# normalBedFile = sys.argv[1]
# pathoBedFile = sys.argv[2]
# outputFile = sys.argv[3]

normalBedFile = snakemake.input[0]
pathoBedFile = snakemake.input[1]
flag = snakemake.params[0]
outputFile = snakemake.output[0]

pathoCNV = defaultdict(list)
with open(pathoBedFile, 'r') as f:
    for x in f:
        x = x.strip().split('\t')
        cnv = "_".join(x[1:4])
        pathoInfo = [x[5], x[6], x[-2]]
        pathoCNV[cnv] = pathoInfo

missing_patho_count = 0

with open(normalBedFile, 'r') as f, open(outputFile, 'w') as g:
    if flag == "before-recurrent":
        print('chrom\tstart\tend\tcn\tcnv\tAS\tDS\tdhfc\tdhbfc\tdhffc\ttools\ttoolNum\tgc\tcnvfilter\tGS\tMS\tNS\tpathogenicity\tPS\tdosageGene', file=g)
    elif flag == "after-recurrent":
        print('chrom\tstart\tend\tcn\tsamples\tsampleNum\taccumScore\tpathogenicity\tpathoScore\tdosageGene', file=g)
    else:
        raise ValueError('You must indicate that this script is used for "before-recurrent" or "after-recurrent".')
    for x in f:
        if x.startswith('chrom'):
            continue
        fileds = x.strip().split('\t')
        cnv = "_".join(fileds[:3])
        if cnv in pathoCNV.keys():
            print(x.strip(), *pathoCNV[cnv], sep='\t', file=g)
        else:
            missing_patho_count += 1
            placeholder = ["Uncertain significance", "0.6", "0.0"]
            print(x.strip(), *placeholder, sep='\t', file=g)

if missing_patho_count > 0:
    import sys
    print(f"WARNING: {missing_patho_count} CNV(s) had no pathogenicity prediction; "
          f"assigned 'Uncertain significance' as default.", file=sys.stderr)
