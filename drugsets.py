# python script to run drug gene set analysis 

# import packages 
import os
import argparse
import subprocess
import numpy as np
import pandas as pd
import drugsets_func as df

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
parser.add_argument('--conditional', '-c', default='no', type=str, choices=['yes','no'],
    help='"yes" will run competitive gene-set analysis in MAGMA while conditioning on a gene set of all druggable genes, "no" will run competitive gene-set analysis without any conditional analysis.',
    required=True)
parser.add_argument('--setsize', '-s', default=2, type=int,
    help='Minimum drug gene set size. Minimum size is 2.',
    required=False)
parser.add_argument('--id', '-i', default='entrez', type=str, choices=['entrez', 'ensembl', 'ensembl92'],
    help='Indicate which gene naming convention is used for your genes.raw file. Options are "entrez" and "ensembl v105", and "ensembl v92". \
        If you ran MAGMA using FUMA, then use "ensembl92"',
    required=False)   
parser.add_argument('--enrich', '-e', default=None, type=str, choices=['atc', 'moa', 'ind', 'all'],
    help='Test drug category for enrichment.',
    required=False)
parser.add_argument('--nsize', '-n', default=5, type=float,
    help = 'Set minimum sample size for drug categories being tested for enrichment.',
    required=False)
parser.add_argument('--showlog', '-l', default='no', type=str, choices = ['no', 'yes'],
    help = 'Print MAGMA output to terminal.',
    required=False)

# parse arguments 
args = parser.parse_args()

# print welcome
print('\n| ----- Welcome to DRUGSETS v1.0 ----- |\n')
print('Reading input...\n')

# check input data
if args.geneassoc[-10:] == '.genes.raw':
    next
else:
    print('ERROR: Gene association file does not end in ".genes.raw". Please check MAGMA gene association input file and try again.')
    quit()

# print input arguments
print('Input arguments used:\n')
for arg in vars(args):
    print('\t', arg,'=', getattr(args, arg))

print('\nRunning drug gene set analysis in MAGMA...\n')

# set base directories
DIR = os.path.dirname(__file__)
DATADIR = os.path.normpath(os.path.join(DIR, 'DATA'))
OUTDIR = os.path.normpath(os.path.join(DIR, 'OUTPUT'))
GENESETDIR = os.path.normpath(os.path.join(DATADIR, 'GENESETS'))
ANNOTDIR = os.path.normpath(os.path.join(DATADIR, 'MAGMA_ANNOT'))

# set filepaths and minimum gene sets size if gene's are named using ENTREZ
if args.id == 'entrez':

    # set gene sets filepaths if setsize is default 
    if args.setsize == 2:
        solo = os.path.normpath(os.path.join(GENESETDIR, 'entrez_genesets.txt'))
        atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_entrez_sets.txt'))
        moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_entrez_sets.txt'))
        ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_entrez_sets.txt'))

    # set file paths for custom minimum gene set size 
    else:
        with open( GENESETDIR+'/entrez_genesets.txt') as oldfile, open( GENESETDIR+'/entrez_genesets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        solo = os.path.normpath(os.path.join(GENESETDIR, 'entrez_genesets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/atc_entrez_sets.txt') as oldfile, open( GENESETDIR+'/atc_entrez_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_entrez_sets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/moa_entrez_sets.txt') as oldfile, open( GENESETDIR+'/moa_entrez_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_entrez_sets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/ind_entrez_sets.txt') as oldfile, open( GENESETDIR+'/ind_entrez_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_entrez_sets_min%d.txt' % args.setsize))

# set filepaths and minimum gene sets size if gene's are named using ENSEMBL
elif args.id == 'ensembl':
    
    # set gene sets filepaths if setsize is default 
    if args.setsize == 2:
        solo = os.path.normpath(os.path.join(GENESETDIR, 'ensembl_genesets.txt'))
        atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_ensembl_sets.txt'))
        moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_ensembl_sets.txt'))
        ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_ensembl_sets.txt'))

    # set file paths for custom minimum gene set size 
    else:
        with open( GENESETDIR+'/ensembl_genesets.txt') as oldfile, open( GENESETDIR+'/ensembl_genesets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        solo = os.path.normpath(os.path.join(GENESETDIR, 'ensembl__genesets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/atc_ensembl_sets.txt') as oldfile, open( GENESETDIR+'/atc_ensembl_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_ensembl_sets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/moa_ensembl_sets.txt') as oldfile, open( GENESETDIR+'/moa_ensembl_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_ensembl_sets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/ind_ensembl_sets.txt') as oldfile, open( GENESETDIR+'/ind_ensembl_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_ensembl_sets_min%d.txt' % args.setsize))

# set filepaths and minimum gene sets size if gene's are named using ENSEMBL
elif args.id == 'ensembl92':
    
    # set gene sets filepaths if setsize is default 
    if args.setsize == 2:
        solo = os.path.normpath(os.path.join(GENESETDIR, 'ensembl_genesets92.txt'))
        atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_ensembl_sets92.txt'))
        moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_ensembl_sets92.txt'))
        ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_ensembl_sets92.txt'))

    # set file paths for custom minimum gene set size 
    else:
        with open( GENESETDIR+'/ensembl_genesets92.txt') as oldfile, open( GENESETDIR+'/ensembl92_genesets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        solo = os.path.normpath(os.path.join(GENESETDIR, 'ensembl92_genesets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/atc_ensembl_sets92.txt') as oldfile, open( GENESETDIR+'/atc_ensembl92_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        atc = os.path.normpath(os.path.join(GENESETDIR, 'atcs_ensembl92_sets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/moa_ensembl_sets92.txt') as oldfile, open( GENESETDIR+'/moa_ensembl92_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_ensembl92_sets_min%d.txt' % args.setsize))

        with open( GENESETDIR+'/ind_ensembl_sets92.txt') as oldfile, open( GENESETDIR+'/ind_ensembl92_sets_min{args.setsize}.txt', 'w') as newfile:
            for line in oldfile:
                if len(line.split('\t')) -3 >= int(args.setsize):
                    newfile.write(line)
        ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_ensembl92_sets_min%d.txt' % args.setsize))


# set MAGMA annotation filepath
annot = os.path.normpath(os.path.join(ANNOTDIR, args.geneassoc))

# set OUTPUT filepath
output = os.path.normpath(os.path.join(OUTDIR, args.out))

# set drug metadata filepath
if args.id == 'entrez':
    metapath = os.path.normpath(os.path.join(DATADIR, 'entrez_meta.pkl'))

elif args.id == 'ensembl':
        metapath = os.path.normpath(os.path.join(DATADIR, 'ensembl105_meta.pkl'))

elif args.id == 'ensembl92':
        metapath = os.path.normpath(os.path.join(DATADIR, 'ensembl92_meta.pkl'))


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
warnings = open(f'{output}.log').read().count('WARNING:')
print('\n\t%s warnings found (see %s.log for details)' % (int(warnings), output))

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

        # ----- ATC CODE ENRICHMENT ----- #

        if args.enrich in ('atc','all'):
            
            # print
            print("ATC enrichment")

            # explode atc
            meta_atc = meta[meta['ATC3'].notna()]
            meta_atc = meta_atc.explode('ATC3')
            
            # test for enrichment of ATC categories 
            resultsATC, resultsBONFatc = df.enrich(meta_atc, gsa_results, 'atc', args.nsize)

            # save all results 
            print('\n\tSaving all results to %s' % (OUTDIR+'/enrich.atc.txt'))            # print 
            resultsPATHatc = os.path.normpath(os.path.join(OUTDIR, 'enrich.atc.txt'))        # set filepath
            resultsATC.to_csv(resultsPATHatc, index = False, sep = '\t')                        # save 

            # add comments
            with open(resultsPATHatc, 'w') as f:
                f.write('# Group = drug category tested for enrichment\n')
                f.write('# MWU = Wilcoxon Mann Whitney U statistic\n')
                f.write('# P = P value from Wilcoxon Mann Whitney U test\n')
                f.write('# AUC = area under the curve\n\n')
            resultsATC.to_csv(resultsPATHatc, index = False, sep = '\t', mode = 'a')

            # save bonferroni corrected results 
            if resultsBONFatc.empty == False:

                #print
                print('\tSaving Bonferroni significant results to %s\n' % (OUTDIR+'/enrich.bonf.atc.txt'))      # print
                bonfPATHatc = os.path.normpath(os.path.join(OUTDIR, 'enrich.bonf.atc.txt'))                        # set filepath
                resultsBONFatc.to_csv(bonfPATHatc, index = False, sep = '\t')                                         # save 

                with open(bonfPATHatc, 'w') as f:
                    f.write('# Group = drug category tested for enrichment\n')
                    f.write('# MWU = Wilcoxon Mann Whitney U statistic\n')
                    f.write('# P = P value from Wilcoxon Mann Whitney U test\n')
                    f.write('# AUC = area under the curve\n\n')
                resultsBONFatc.to_csv(bonfPATHatc, index = False, sep = '\t', mode = 'a')
        

        # ---- CLINICAL INDICATION ENRICHMENT ----- #         

        if args.enrich in ('ind','all'):

            # print
            print('Clinical indication enrichment')

            # explode indication
            meta_ind = meta.explode('indication')

            # test for enrichment of ATC categories 
            resultsIND, resultsBONFind = df.enrich(meta_ind, gsa_results, 'ind', args.nsize)

            # save all results 
            print('\n\tSaving all results to %s' % (OUTDIR+'/enrich.ind.txt'))            # print 
            resultsPATHind = os.path.normpath(os.path.join(OUTDIR, 'enrich.ind.txt'))        # set filepath
            resultsIND.to_csv(resultsPATHind, index = False, sep = '\t')                        # save 

            # add comments
            with open(resultsPATHind, 'w') as f:
                f.write('# Group = drug category tested for enrichment\n')
                f.write('# MWU = Wilcoxon Mann Whitney U statistic\n')
                f.write('# P = P value from Wilcoxon Mann Whitney U test\n')
                f.write('# AUC = area under the curve\n\n')
            resultsIND.to_csv(resultsPATHind, index = False, sep = '\t', mode = 'a')

            # save bonferroni corrected results 
            if resultsBONFind.empty == False:

                #print
                print('\tSaving Bonferroni significant results to %s\n' % (OUTDIR+'/enrich.bonf.ind.txt'))                 # print
                bonfPATHind = os.path.normpath(os.path.join(OUTDIR, 'enrich.bonf.ind.txt'))                                   # set filepath
                resultsBONFind.to_csv(bonfPATHind, index = False, sep = '\t')                                                    # save 

                # add comments 
                with open(bonfPATHind, 'w') as f:
                    f.write('# Group = drug category tested for enrichment\n')
                    f.write('# MWU = Wilcoxon Mann Whitney U statistic\n')
                    f.write('# P = P value from Wilcoxon Mann Whitney U test\n')
                    f.write('# AUC = area under the curve\n\n')
                resultsBONFind.to_csv(bonfPATHind, index = False, sep = '\t', mode = 'a')

        # ---- MOA ENRICHMENT ----- #         

        if args.enrich in ('moa','all'):

            # print
            print("Mechanism of action enrichment")

            # explode moa
            meta_moa = meta.explode('moa')

            # test for enrichment of ATC categories 
            resultsMOA, resultsBONFmoa = df.enrich(meta_moa, gsa_results, 'moa', args.nsize)
            
            # save all results 
            print('\n\tSaving all results to %s' % (OUTDIR+'/enrich.moa.txt'))                        # print 
            resultsPATHmoa = os.path.normpath(os.path.join(OUTDIR, 'enrich.moa.txt'))                    # set filepath
            resultsMOA.to_csv(resultsPATHmoa, index = False, sep = '\t')                                    # save 

            # add comments
            with open(resultsPATHmoa, 'w') as f:
                f.write('# Group = drug category tested for enrichment\n')
                f.write('# MWU = Wilcoxon Mann Whitney U statistic\n')
                f.write('# P = P value from Wilcoxon Mann Whitney U test\n')
                f.write('# AUC = area under the curve\n\n')
            resultsMOA.to_csv(resultsPATHmoa, index = False, sep = '\t', mode = 'a')

            # save bonferroni corrected results 
            if resultsBONFmoa.empty == False:

                #print
                print('\tSaving Bonferroni significant results to %s\n' % (OUTDIR+'/enrich.bonf.moa.txt'))                 # print
                bonfPATHmoa = os.path.normpath(os.path.join(OUTDIR, 'enrich.bonf.moa.txt'))                                   # set filepath
                resultsBONFmoa.to_csv(bonfPATHmoa, index = False, sep = '\t')                                                    # save 

                # add comments 
                with open(bonfPATHmoa, 'w') as f:
                    f.write('# Group = drug category tested for enrichment\n')
                    f.write('# MWU = Wilcoxon Mann Whitney U statistic\n')
                    f.write('# P = P value from Wilcoxon Mann Whitney U test\n')
                    f.write('# AUC = area under the curve\n\n')
                resultsBONFmoa.to_csv(bonfPATHmoa, index = False, sep = '\t', mode = 'a')

        # remove new gene set files if created new 
        # if args.setsize == 2:
        #     next
        # else:
        #     subprocess.run('rm %s/*min%d.txt' % (GENESETDIR, args.setsize), shell=True)

        # print finished 
        print('\n\tEnrichment analysis finished.\n')

    else: 
        print('To test for enrichment "-drugsets" must be set to "solo".')


            