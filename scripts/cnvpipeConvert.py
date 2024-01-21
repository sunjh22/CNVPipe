#! /usr/bin/env python

# Give overall score for each CNV based on previous metrics

import sys
import pandas as pd
import kmeans1d

# inputFile = sys.argv[1]
# outputFile = sys.argv[2]

inputFile = snakemake.input[0]
outputFile = snakemake.output[0]

cnvList = pd.read_csv(inputFile, sep='\t')

# assign a score based on the predicted pathogenicity of CNV
patho_score = {'Pathogenic':1, 'Likely pathogenic':0.8, 'Uncertain significance':0.6, 'Likely benign': 0.4, 'Benign':0.2}
cnvList['patho_score'] = [patho_score[x] for x in cnvList['pathogenicity']]

# a large dosage score will be assigned if there is a dosage gene in the CNV
cnvList['dosage_score'] = [1 if x else 0 for x in cnvList['dosageGene']]

# scale the other four features: too_num, accumu_score, duphold_score, good_score
cnvList_part = cnvList.loc[:, ['toolNum', 'AS', 'DS', 'GS']]
cnvList_part_scale = (cnvList_part - cnvList_part.min()) / (cnvList_part.max() - cnvList_part.min())
cnvList_part_scale['GS'].fillna(1, inplace=True)
cnvList_part_scale.rename(columns={'toolNum':'toolNum_scale', 'AS':'AS_scale', 'DS':'DS_scale', 'GS':'GS_scale'}, inplace=True)

# construct a new dataframe with scaled features
new_cnvList = pd.concat([cnvList, cnvList_part_scale], axis=1)

# calculate overall score, 'patho_score' has large weight
overall_score = new_cnvList['toolNum_scale'] + new_cnvList['AS_scale'] + new_cnvList['DS_scale'] + new_cnvList['GS_scale'] + new_cnvList['patho_score']*7 + new_cnvList['dosage_score']

# use one-dimensional clustering to separate CNVs into five categories
clusters, centroids = kmeans1d.cluster(overall_score, 5)

# assign a label (star) for each cluster
label_dict = {'0':'*','1':'**','2':'***','3':'****','4':'*****'}
priority_label = [label_dict[str(x)] for x in clusters]
new_cnvList.insert(15, 'overall_score', overall_score.round(3))

# assign labels by setting specific quantile cutoff
# cnv_quantiles = new_cnvList['overall_score'].quantile([0.5,0.75])
# priority_label = []
# for x in new_cnvList['overall_score']:
#     if x <= cnv_quantiles[0.5]:
#         priority_label.append('*')
#     elif cnv_quantiles[0.5] < x <= cnv_quantiles[0.75]:
#         priority_label.append('**')
#     else:
#         priority_label.append('***')

new_cnvList.insert(16, 'priority', priority_label)
new_cnvList.sort_values(by='overall_score', axis=0, ascending=False, inplace=True)
# print(new_cnvList.loc[:,['chrom','start','end','cn','cnv','toolNum','pathogenicity','dosageGene','overall_score','priority']].head())
# new_cnvList.loc[:,'chrom':'priority'].to_csv(outputFile, sep='\t', index=False)
new_cnvList.loc[:,['chrom','start','end','cn','cnv','toolNum','pathogenicity','dosageGene','overall_score','priority']].to_csv(outputFile, sep='\t', index=False)