#! /usr/bin/env python

# Give overall score for each CNV based on previous metrics

import sys
import pandas as pd

# inputFile = sys.argv[1]
# outputFile = sys.argv[2]

inputFile = snakemake.input[0]
outputFile = snakemake.output[0]

cnvList = pd.read_csv(inputFile, sep='\t')
patho_score = {'Pathogenic':1, 'Likely pathogenic':0.8, 'Uncertain significance':0.6, 'Likely benign': 0.4, 'Benign':0.2}
cnvList['patho_score'] = [patho_score[x] for x in cnvList['pathogenicity']]
cnvList['cnvfilter_score'] = [1 if x else 0 for x in cnvList['cnvfilter']]
cnvList['dosage_score'] = [1 if x else 0 for x in cnvList['dosageGene']]
cnvList_part = cnvList.loc[:, ['toolNum', 'AS', 'DS', 'GS']]
cnvList_part_scale = (cnvList_part - cnvList_part.min()) / (cnvList_part.max() - cnvList_part.min())
cnvList_part_scale.rename(columns={'toolNum':'toolNum_scale', 'AS':'AS_scale', 'DS':'DS_scale', 'GS':'GS_scale'}, inplace=True)
new_cnvList = pd.concat([cnvList, cnvList_part_scale], axis=1)
overall_score = new_cnvList['toolNum_scale']*2 + new_cnvList['AS_scale'] + new_cnvList['DS_scale'] + new_cnvList['GS_scale'] + new_cnvList['cnvfilter_score'] + new_cnvList['patho_score'] + new_cnvList['dosage_score']
new_cnvList.insert(15, 'overall_score', overall_score)
cnv_quantiles = new_cnvList['overall_score'].quantile([0.5,0.75])
priority_label = []
for x in new_cnvList['overall_score']:
    if x <= cnv_quantiles[0.5]:
        priority_label.append('*')
    elif cnv_quantiles[0.5] < x <= cnv_quantiles[0.75]:
        priority_label.append('**')
    else:
        priority_label.append('***')
new_cnvList.insert(16, 'priority_label', priority_label)
new_cnvList.sort_values(by='overall_score', axis=0, ascending=False, inplace=True)
# print(new_cnvList.loc[:,'chrom':'priority_label'].head())
new_cnvList.loc[:,'chrom':'priority_label'].to_csv(outputFile, sep='\t', index=False)