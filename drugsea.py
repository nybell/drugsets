# python script to run drug gene set analysis 

# import packages 
import os
import argparse
import subprocess
import numpy as np
import pandas as pd
import drugsea_func as df

# parse arguments 
parser = argparse.ArgumentParser()
parser.add_argument('--geneassoc', '-g', default=None, type=str,
    help='Filename of gene associations from MAGMA (.genes.raw).',
    required=True)
parser.add_argument('--drugsets', '-d', default='solo', type=str, choices=['solo', 'atc', 'moa', 'ind'],
    help='Type of drug gene set to use (individual, ATC code, mechanism of action, clinical indication).',
    required=True)
parser.add_argument('--out', '-o', default=None, type=str,
    help='Filename of output.',
    required=True)
parser.add_argument('--setsize', '-s', default=2, type=int,
    help='Minimum drug gene set size. Minimum size is 2.',
    required=False)
parser.add_argument('--enrich', '-e', default=None, type=str, choices=['atc', 'moa', 'ind'],
    help='Test drug category for enrichment.',
    required=False)
parser.add_argument('--nsize', '-n', default=5, type=float,
    help = 'Set minimum sample size for drug categories being tested for enrichment.',
    required=False)
parser.add_argument('--showlog', '-l', default='no', type=str, choices = ['no', 'yes'],
    help = 'Print MAGMA output to terminal.',
    required=False)

# print 
print('\n| ----- Welcome to DRUGSEA v1.0 ----- |\n')
print('Reading input...\n')
print('Running drug gene set analysis in MAGMA...\n')

# parse arguments
args = parser.parse_args()

# set base directories
DIR = os.path.dirname(__file__)
DATADIR = os.path.normpath(os.path.join(DIR, 'DATA'))
OUTDIR = os.path.normpath(os.path.join(DIR, 'OUTPUT'))
GENESETDIR = os.path.normpath(os.path.join(DATADIR, 'GENESETS'))
ANNOTDIR = os.path.normpath(os.path.join(DATADIR, 'MAGMA_ANNOT'))

# set gene sets filepaths if setsize is default 
if args.setsize == 2:
    solo = os.path.normpath(os.path.join(GENESETDIR, 'drug_genesets.txt'))
    atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_sets.txt'))
    moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_sets.txt'))
    ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_sets.txt'))

# set file paths for custom minimum gene set size 
else:
    subprocess.run("awk 'NF>=%d' %s > %s" % (args.setsize+1, GENESETDIR+'/drug_genesets.txt', GENESETDIR+'/drugsets_min%d.txt' % args.setsize), shell=True)
    solo = os.path.normpath(os.path.join(GENESETDIR, 'drugsets_min%d.txt' % args.setsize))

    subprocess.run('awk "NF>=%d" %s > %s' % (args.setsize+1, GENESETDIR+'/atc_sets.txt', GENESETDIR+'/atcsets_min%d.txt' % args.setsize), shell=True)
    atc = os.path.normpath(os.path.join(GENESETDIR, 'atcsets_min%d.txt' % args.setsize))

    subprocess.run('awk "NF>=%d" %s > %s' % (args.setsize+1, GENESETDIR+'/moa_sets.txt', GENESETDIR+'/moasets_min%d.txt' % args.setsize), shell=True)
    moa = os.path.normpath(os.path.join(GENESETDIR, 'moasets_min%d.txt' % args.setsize))

    subprocess.run('awk "NF>=%d" %s > %s' % (args.setsize+1, GENESETDIR+'/ind_sets.txt', GENESETDIR+'/indsets_min%d.txt' % args.setsize), shell=True)
    ind = os.path.normpath(os.path.join(GENESETDIR, 'indsets_min%d.txt' % args.setsize))


# set MAGMA annotation filepath
annot = os.path.normpath(os.path.join(ANNOTDIR, args.geneassoc))

# set OUTPUT filepath
output = os.path.normpath(os.path.join(OUTDIR, args.out))

# set drug metadata filepath
metapath = os.path.normpath(os.path.join(DATADIR, 'geneset_meta.pkl'))

# execute gene set analysis 
if args.drugsets == 'solo':
    if args.showlog == 'no':
        df.run_task('magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, solo, output))
    else:
        subprocess.run(['magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, solo, output)], \
            shell = True, check = True)

elif args.drugsets == 'atc':
    if args.showlog == 'no':
        df.run_task('magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, atc, output))
    else:
        subprocess.run(['magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, atc, output)], \
            shell = True, check = True)

elif args.drugsets == 'moa':
    if args.showlog == 'no':
        df.run_task('magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, moa, output))
    else:
        subprocess.run(['magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, moa, output)], \
            shell = True, check = True)

elif args.drugsets == 'ind':
    if args.showlog == 'no':
        df.run_task('magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, ind, output))
    else:
        subprocess.run(['magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, ind, output)], \
            shell = True, check = True)

# print log 
warnings = subprocess.run('grep WARNING: ./OUTPUT/%s.log | wc -l ' % args.out, shell = True, capture_output=True)
print('\n\t%s warnings found (see %s.log for details)' % (int(warnings.stdout), output))

# print result locations 
print('\tResults for all drug gene sets saving to %s' % (OUTDIR+'/%s.gsa.out' % args.out))
print('\tResults for significant drug gene sets saving to %s\n' % (OUTDIR+'/%s.gsa.set.genes.out' % args.out))


# print done
print('\tDrug gene set analysis finished.\n')


# enrichment analysis
if args.enrich is not None:

    if args.drugsets == 'solo':
        
        #print 
        print('Running enrichment analysis...\n\n')

        # set file path for .gsa.out results file
        gsa = (output+'.gsa.out')

        # load gsa results 
        gsa_results = pd.read_csv(gsa, delimiter= "\s+", comment='#') 

        # load drug meta data
        meta = pd.read_pickle(metapath)              

        # explode moa
        meta_moa = meta.explode('moa')

        # ----- ATC CODE ENRICHMENT ----- #

        if args.enrich == 'atc':

            # explode atc
            meta_atc = meta[meta['Therapeutic_classification_level_III_ATC'].notna()]
            meta_atc = meta_atc.explode('Therapeutic_classification_level_III_ATC')
            
            # test for enrichment of ATC categories 
            results, resultsBONF = df.enrich(meta_atc, gsa_results, 'atc', args.nsize)

            # save all results 
            print('\n\tSaving all results to %s' % (OUTDIR+'/enrich.%s.txt' % args.enrich))            # print 
            resultsPATH = os.path.normpath(os.path.join(OUTDIR, 'enrich.%s.txt' % args.enrich))        # set filepath
            results.to_csv(resultsPATH, index = False, sep = '\t')                                     # save 

            # save bonferroni corrected results 
            if resultsBONF.empty == False:

                #print
                print('\tSaving Bonferroni significant results to %s' % (OUTDIR+'/enrich.bonf.%s.txt' % args.enrich))      # print
                bonfPATH = os.path.normpath(os.path.join(OUTDIR, 'enrich.bonf.%s.txt' % args.enrich))                      # set filepath
                resultsBONF.to_csv(bonfPATH, index = False, sep = '\t')                                                    # save 


        # ---- CLINICAL INDICATION ENRICHMENT ----- #         

        elif args.enrich == 'ind':

            # explode indication
            meta_ind = meta.explode('indication')

            # test for enrichment of ATC categories 
            results, resultsBONF = df.enrich(meta_ind, gsa_results, 'ind', args.nsize)

            # save all results 
            print('\n\tSaving all results to %s' % (OUTDIR+'/enrich.%s.txt' % args.enrich))            # print 
            resultsPATH = os.path.normpath(os.path.join(OUTDIR, 'enrich.%s.txt' % args.enrich))        # set filepath
            results.to_csv(resultsPATH, index = False, sep = '\t')                                     # save 

            # save bonferroni corrected results 
            if resultsBONF.empty == False:

                #print
                print('\tSaving Bonferroni significant results to %s' % (OUTDIR+'/enrich.bonf.%s.txt' % args.enrich))      # print
                bonfPATH = os.path.normpath(os.path.join(OUTDIR, 'enrich.bonf.%s.txt' % args.enrich))                      # set filepath
                resultsBONF.to_csv(bonfPATH, index = False, sep = '\t')                                                    # save 

        # ---- MOA ENRICHMENT ----- #         

        elif args.enrich == 'moa':

            # explode moa
            meta_moa = meta.explode('moa')

            # test for enrichment of ATC categories 
            results, resultsBONF = df.enrich(meta_moa, gsa_results, 'moa', args.nsize)

            # save all results 
            print('\n\tSaving all results to %s' % (OUTDIR+'/enrich.%s.txt' % args.enrich))                        # print 
            resultsPATH = os.path.normpath(os.path.join(OUTDIR, 'enrich.%s.txt' % args.enrich))                    # set filepath
            results.to_csv(resultsPATH, index = False, sep = '\t')                                                 # save 

            # save bonferroni corrected results 
            if resultsBONF.empty == False:

                #print
                print('\tSaving Bonferroni significant results to %s' % (OUTDIR+'/enrich.bonf.%s.txt' % args.enrich))                   # print
                bonfPATH = os.path.normpath(os.path.join(OUTDIR, 'enrich.bonf.%s.txt' % args.enrich))                                   # set filepath
                resultsBONF.to_csv(bonfPATH, index = False, sep = '\t')                                                                 # save 

        # remove new gene set files if created new 
        if args.setsize == 2:
            next
        else:
            subprocess.run('rm %s/*min%d.txt' % (GENESETDIR, args.setsize), shell=True)

        # print finished 
        print('\n\tEnrichment analysis finished.\n')

    else: 
        print('To test for enrichment "-drugsets" must be set to "solo".')


            
