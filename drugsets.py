# python script to run drug gene set analysis 

# -------------------------------------------------------------------------- #
##### ----- PART 1: IMPORT PACKAGES, PARSE ARGUMENTS, CHECK INPUTS ----- #####
# -------------------------------------------------------------------------- #

# import packages 
import os
import argparse
import pandas as pd
import drugsets_func as df

# parse arguments 
parser = argparse.ArgumentParser()
parser.add_argument('--geneassoc', '-g', default=None, type=str,
    help='Filename of gene associations from MAGMA (.genes.raw).',
    required=True)
parser.add_argument('--drugsets', '-d', default='solo', type=str, choices=['solo', 'atc', 'moa', 'ind', 'all'],
    help='Type of drug gene set to use (individual, ATC code, mechanism of action, clinical indication).',
    required=True)
parser.add_argument('--out', '-o', default=None, type=str,
    help='File path + prefix name of output file.',
    required=True)
parser.add_argument('--conditional', '-c', default='yes', type=str, choices=['yes','no'],
    help='"yes" will run competitive gene-set analysis in MAGMA while conditioning on a gene set of all druggable genes, "no" will run competitive gene-set analysis without any conditional analysis',
    required=False)
parser.add_argument('--setsize', '-s', default=2, type=int,
    help='Minimum drug gene set size. Minimum size is 2.',
    required=False)
parser.add_argument('--id', '-i', default='entrez', type=str, choices=['entrez', 'ensembl', 'ensembl92'],
    help='Indicate which gene naming convention is used for your genes.raw file. Options are "entrez" and "ensembl v105", and "ensembl v92". \
        If you ran MAGMA using FUMA, then use "ensembl92"',
    required=False)   
parser.add_argument('--enrich', '-e', default=None, type=str, choices=['atc', 'moa', 'ind', 'all'],
    help='Test drug category for enrichment using linear regression model.',
    required=False)
parser.add_argument('--enrich_auc', '-auc', default=None, type=str, choices=['atc', 'moa', 'ind', 'all'],
    help='Test drug category for enrichment using drug AUC.',
    required=False)
parser.add_argument('--nsize', '-n', default=5, type=float,
    help = 'Set minimum sample size for drug categories being tested for enrichment.',
    required=False)
parser.add_argument('--correct', '-p', default='bonf', type=str, choices=['bonf', 'fdr'],
    help = 'Select correction to be used for multiple testing for drug group enrichment: Bonferroni or FDR',
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

# ------------------------------------------------------------------------------------------------- #
##### ----- PART 2: DEFINE FILEPATHS AND GENESETS BASED ON ID, SIZE, AND CONDITION INPUTS ----- #####
# ------------------------------------------------------------------------------------------------- #

# output filepath
out_head_tail = os.path.split(args.out)
out_path = out_head_tail[0]              # define output file path 
out_name = out_head_tail[1]              # define output file name

# check out path
if out_path == '':
    out_path = os.getcwd()         # if no file path is given for the output it will put it all in the current working directory
    output = out_path+'/'+out_name
elif '.' in out_path:
    out_path = out_path.replace('.', os.getcwd())       # replace '.' with full file path to avoid errors later on
    output = out_path+'/'+out_name
else:
    output = out_path+'/'+out_name

# Mapping conditions to prefixes
cond_prefixes = {'yes': '_cond', 'no': ''}

# Setting directories
DIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.normpath(os.path.join(DIR, 'DATA'))
OUTDIR = out_path
GENESETDIR = os.path.normpath(os.path.join(DATADIR, 'GENESETS'))
ANNOTDIR = os.path.normpath(os.path.join(DATADIR, 'MAGMA_ANNOT'))

# Filepath setup
gene_id_conditions = ['entrez', 'ensembl', 'ensembl92']

filepaths = {}

for id in gene_id_conditions:
    filepaths[id] = {
        'solo': df.set_paths(f"solo_{id}", f"{cond_prefixes[args.conditional]}", args.setsize, GENESETDIR),
        'atc': df.set_paths(f"atc_{id}", f"{cond_prefixes[args.conditional]}", args.setsize, GENESETDIR),
        'moa': df.set_paths(f"moa_{id}", f"{cond_prefixes[args.conditional]}", args.setsize, GENESETDIR),
        'ind': df.set_paths(f"ind_{id}", f"{cond_prefixes[args.conditional]}", args.setsize, GENESETDIR)
    }

solo = filepaths[args.id]['solo']
atc = filepaths[args.id]['atc']
moa = filepaths[args.id]['moa']
ind = filepaths[args.id]['ind']

# set MAGMA annotation filepath
if '\\' in args.geneassoc:
    print(not'/' in args.geneassoc or not '\\' in args.geneassoc)

if './' in args.geneassoc or '.\\' in args.geneassoc:
    annot = args.geneassoc.replace('.', os.getcwd(),1)
elif '/' in args.geneassoc or '\\' in args.geneassoc:                        # if no file path is given for the .genes.raw file it assumes it is in the working directory
    annot = args.geneassoc
else:
    annot = os.path.normpath(os.getcwd() + '/' + args.geneassoc)

# set path for meta data pkl files
if args.id == 'entrez':
    metapath_pkl = os.path.normpath(os.path.join(DATADIR, 'entrez_meta.pkl'))
elif args.id == 'ensembl':
    metapath_pkl = os.path.normpath(os.path.join(DATADIR, 'ensembl105_meta.pkl'))
elif args.id == 'ensembl92':
    metapath_pkl = os.path.normpath(os.path.join(DATADIR, 'ensembl92_meta.pkl'))

    
# ------------------------------------------------------ #
##### ----- PART 3: RUN DRUG GENE SET ANALYSIS ----- #####
# ------------------------------------------------------ #

# individual drug gene set analysis
#Set output labels for each chosen analysis
if args.drugsets == 'all':
    analysis_labels = ['_SOLO','_ATC', '_MOA', '_IND']
else:
    analysis_labels = ['_' + args.drugsets.upper()]

#Loop over each analysis output label and match input file for corresponding analysis
for label  in analysis_labels:
    if label == '_SOLO':
        analysis = solo
    elif label == '_ATC':
        analysis = atc
    elif label == '_MOA':
        analysis = moa
    elif label == '_IND':
        analysis = ind

    if args.conditional == 'no':
        print('\nRunning %s drug gene set analysis in MAGMA...\n' % (label.strip('_')))
        df.run_task('magma --gene-results %s --set-annot %s --settings gene-info --out %s' % (annot, analysis, os.path.normpath((output + label))))
        
        # print log 
        warnings = open(os.path.normpath(f'{output+label}.log')).read().count('WARNING:')
        print('\n\t%s warnings found (see %s.log for details)' % (int(warnings), os.path.normpath((output + label))))

        # print result locations 
        print('\tResults for all drug gene sets saving to %s' % (os.path.normpath(OUTDIR+'/%s.gsa.out' % (out_name + label))))
        print('\tResults for significant drug gene sets saving to %s\n' % (os.path.normpath(OUTDIR+'/%s.gsa.set.genes.out' % (out_name + label))))

    elif args.conditional == 'yes':
        
        print('\nRunning conditional %s drug gene set analysis in MAGMA...\n' % (label.strip('_')))
        df.run_task('magma --gene-results %s --set-annot %s --settings gene-info --model condition=druggable --out %s' % (annot, analysis, os.path.normpath((output + label))))
        
        # print log 
        warnings = open(os.path.normpath(f'{output+label}.log')).read().count('WARNING:')
        print('\n\t%s warnings found (see %s.log for details)' % (int(warnings), os.path.normpath((output + label))))

        # print result locations 
        print('\tResults for all drug gene sets saving to %s' % (os.path.normpath(OUTDIR+'/%s.gsa.out' % (out_name + label))))
        print('\tResults for significant drug gene sets saving to %s\n' % (os.path.normpath(OUTDIR+'/%s.gsa.set.genes.out' % (out_name + label))))

# print done
print('Drug gene set analysis finished.\n')

# -------------------------------------------------------------- #
##### ----- PART 4: DRUG GROUP ENRICHMENT - REGRESSION ----- #####
# -------------------------------------------------------------- #

# enrichment analysis
if args.enrich is not None:

    if (args.drugsets == 'solo' or args.drugsets == 'all'):

        # set full file paths for .raw file, gene set file
        working = os.path.dirname(os.path.abspath(__file__)) + '/'
        
        #print 
        print('Running %s enrichment analysis...\n\n' % (args.enrich.upper()))

        # set file path for .gsa.out results file
        if os.path.exists(os.path.normpath((output+'_SOLO'+'.gsa.out'))):
            gsa = os.path.normpath((output+'_SOLO'+'.gsa.out'))
        elif os.path.exists(os.path.normpath((output+'_SOLO'+'.gsa.out.txt'))):
            gsa = os.path.normpath((output+'_SOLO'+'.gsa.out.txt'))              

        # compute covariance 
        print('\tComputing correlation matrix...')
        df.run_task_r('Rscript --vanilla %s %s %s %s %s %s' % (os.path.normpath(working+'compute_corrs.R'), annot, os.path.normpath(solo), (gsa), os.path.normpath('/'+out_name), OUTDIR)) 

        # define filepath to set.corrs.rdata and to metadata.rdata file
        corrdata = os.path.normpath(out_path+'/'+'%s_setcorrs.rdata' % out_name)
        metaRdata = os.path.normpath(working+'DATA/metadata.rdata')

        # run either one group or all 

        if args.enrich != 'all':
            
            # compute dependent linear regression
            print('\tRunning dependent linear regression model...')
            df.run_task_r('Rscript --vanilla %s %s %s %s %s %s %s %s' % (os.path.normpath(working+'compute_lnreg.R'), corrdata, metaRdata, args.enrich.lower(), args.nsize, out_name, (OUTDIR+os.path.sep), args.correct))

        elif args.enrich == 'all':
            
            # define groups to loop through 
            groups = ['atc','moa','ind']

            # loop though types of drug groups and run dependent linear regression model for each type
            for g in groups:

                # compute dependent linear regression
                print('\tRunning dependent linear regression model for %s groups...' % g.upper())
                df.run_task_r('Rscript --vanilla %s %s %s %s %s %s %s %s' % (os.path.normpath(working+'compute_lnreg.R'), corrdata, metaRdata, g, args.nsize, (out_name+'_'+g.upper()), (OUTDIR+os.path.sep), args.correct))

        # remove correlation matrix file
        os.remove(corrdata)

        # remove new gene set files if created new 
        # if args.setsize == 2:
        #     next
        # else:
        #     subprocess.run('rm %s/*min%d.txt' % (GENESETDIR, args.setsize), shell=True)

        # print finished 
        print('\nEnrichment analysis finished.\n')

    else: 
        print('To test for enrichment "-drugsets" must be set to "solo".')


# ------------------------------------------------------- #
##### ----- PART 5: DRUG GROUP ENRICHMENT - AUC ----- #####
# ------------------------------------------------------- #

if args.enrich_auc is not None:
    if args.drugsets == 'solo':
        
        #print 
        print('Running %s drug AUC analysis...\n\n' % (args.enrich_auc.upper()))
        
        # set file path for .gsa.out results file
        gsa = (output+'_SOLO.gsa.out')

        # load gsa results 
        gsa_results = pd.read_csv(gsa, delimiter= "\s+", comment='#') 
        
        # load drug meta data
        meta = pd.read_pickle(metapath_pkl)              
        
        # ----- ATC CODE ENRICHMENT ----- #
        if args.enrich_auc in ('atc','all'):
            
            # print
            print("ATC enrichment")
            # explode atc
            meta_atc = meta[meta['ATC3'].notna()]
            meta_atc = meta_atc.explode('ATC3')
            
            # test for enrichment of ATC categories 
            resultsATC, resultsBONFatc = df.enrich_auc(meta_atc, gsa_results, 'atc', args.nsize)

            # save all results 
            resultsPATHatc = os.path.normpath(os.path.join(OUTDIR, out_name+'.enrich_auc.atc.txt'))        # set filepath
            print('\n\tSaving all results to %s' % (resultsPATHatc))                                       # print 
            resultsATC.to_csv(resultsPATHatc, index = False, sep = '\t')                                   # save 

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
                bonfPATHatc = os.path.normpath(os.path.join(OUTDIR, out_name+'.enrich_auc.bonf.atc.txt'))        # set filepath
                print('\tSaving Bonferroni significant results to %s\n' % (bonfPATHatc))                          # print
                resultsBONFatc.to_csv(bonfPATHatc, index = False, sep = '\t')                                     # save 

                # add comments
                with open(bonfPATHatc, 'w') as f:
                    f.write('# Group = drug category tested for enrichment\n')
                    f.write('# MWU = Wilcoxon Mann Whitney U statistic\n')
                    f.write('# P = P value from Wilcoxon Mann Whitney U test\n')
                    f.write('# AUC = area under the curve\n\n')
                resultsBONFatc.to_csv(bonfPATHatc, index = False, sep = '\t', mode = 'a')
        
        # ---- CLINICAL INDICATION ENRICHMENT ----- #         
        if args.enrich_auc in ('ind','all'):

            # print
            print('Clinical indication enrichment')

            # explode indication
            meta_ind = meta.explode('indication')

            # test for enrichment of ATC categories 
            resultsIND, resultsBONFind = df.enrich_auc(meta_ind, gsa_results, 'ind', args.nsize)

            # save all results 
            resultsPATHind = os.path.normpath(os.path.join(OUTDIR, out_name+'.enrich_auc.ind.txt'))        # set filepath
            print('\n\tSaving all results to %s' % (resultsPATHind))                          # print 
            resultsIND.to_csv(resultsPATHind, index = False, sep = '\t')                                    # save 

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
                bonfPATHind = os.path.normpath(os.path.join(OUTDIR, out_name+'.enrich_auc.bonf.ind.txt'))                       # set filepath
                print('\tSaving Bonferroni significant results to %s\n' % (bonfPATHind))        # print
                resultsBONFind.to_csv(bonfPATHind, index = False, sep = '\t')                                         # save 

                # add comments 
                with open(bonfPATHind, 'w') as f:
                    f.write('# Group = drug category tested for enrichment\n')
                    f.write('# MWU = Wilcoxon Mann Whitney U statistic\n')
                    f.write('# P = P value from Wilcoxon Mann Whitney U test\n')
                    f.write('# AUC = area under the curve\n\n')
                resultsBONFind.to_csv(bonfPATHind, index = False, sep = '\t', mode = 'a')

        # ---- MOA ENRICHMENT ----- #         
        if args.enrich_auc in ('moa','all'):

            # print
            print("Mechanism of action enrichment")

            # explode moa
            meta_moa = meta.explode('moa')

            # test for enrichment of ATC categories 
            resultsMOA, resultsBONFmoa = df.enrich_auc(meta_moa, gsa_results, 'moa', args.nsize)
            
            # save all results 
            resultsPATHmoa = os.path.normpath(os.path.join(OUTDIR, out_name+'.enrich_auc.moa.txt'))                 # set filepath
            print('\n\tSaving all results to %s' % (resultsPATHmoa))                        # print 
            resultsMOA.to_csv(resultsPATHmoa, index = False, sep = '\t')                              # save 

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
                bonfPATHmoa = os.path.normpath(os.path.join(OUTDIR, out_name+'.enrich_auc.bonf.moa.txt'))                                # set filepath
                print('\tSaving Bonferroni significant results to %s\n' % (bonfPATHmoa))                 # print
                resultsBONFmoa.to_csv(bonfPATHmoa, index = False, sep = '\t')                                              # save 

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
            