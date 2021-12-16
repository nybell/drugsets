#drugsea_func.py>

# import packages 
import os
import sys
import argparse
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import stats
# import matplotlib.pyplot as plt 
from sklearn.metrics import roc_curve, auc
from subprocess import Popen, PIPE, CalledProcessError

# function to test for enrichment of significant drug using Wilcoxon Mann Whitney U test
def testMWU(data, group, subgroup):     # inputs: data = data; set subgroup = name of subgroup (e.g., 'N05A' within all ATC groups); group = string of group (e.g., 'ATC3')
    x = data[data[group]==subgroup]
    x = x['P']
    y = data[data[group]!=subgroup]
    y = y['P']
    return stats.mannwhitneyu(x, y, use_continuity=True, alternative='greater', axis=0, method='auto')

def drugAUC(data, group, subgroup):
    y_score = data['P']
    y_true = data[group]==subgroup                                      
    y_true = y_true*1
    fpr, tpr, threshold = roc_curve(y_true,y_score)
    roc_auc = auc(fpr,tpr)
    return roc_auc

def enrich(data, results, group, nsize):
            # merge MAGMA results data with drug meta data 
            mergedData = data.merge(results, how = "right", left_on = "DRUG", right_on = "VARIABLE")
            mergedData['P'] = abs(np.log10(mergedData['P']))
            mergedData = mergedData.sort_values(by=['P'], ascending=False)

            # drop rows with NA drug names and explode ATC III classifications 
            mergedData = mergedData[mergedData['DRUG'].notna()]

            if group == 'atc':
                mergedData = mergedData.explode('Therapeutic_classification_level_III_ATC')
                mergedData = mergedData.rename(columns={'Therapeutic_classification_level_III_ATC':'atc'})

            if group == 'ind':
                mergedData = mergedData.rename(columns={'indication':'ind'})


            # ATC codes 
            codes = pd.DataFrame(mergedData[group].value_counts())
            codes.reset_index(level=0, inplace=True)
            codes = codes[codes[group] >= nsize]

            # run testMWU for every ATC III code with more than 5 drugs and then store in variable 
            results = []

            for code in tqdm(codes['index']):

                # test for enrichment using WMW U test
                mwu, p = testMWU(mergedData, group, code)      

                # calculate AUC enrichment curves
                roc_auc = drugAUC(mergedData, group, code)

                # append to results
                result = [code,mwu,p, roc_auc]
                results.append(result)

            # convert to csv file
            results = pd.DataFrame(results, columns = ['GROUP', 'MWU','P', 'AUC'])

            # bonferroni correct ATC results 
            bonf = 0.05/len(codes)
            resultsBONF = results[results['P'] < bonf]

            # return results 
            return results, resultsBONF

def run_task(cmd):
    with Popen(cmd, shell=True, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for b in tqdm(p.stdout):
            # print(b, end='') # b is the byte from stdout
            next

    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

