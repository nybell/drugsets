# python script to run drug gene set analysis 

# -------------------------------------------------------------------------- #
##### ----- PART 1: IMPORT PACKAGES, PARSE ARGUMENTS, CHECK INPUTS ----- #####
# -------------------------------------------------------------------------- #

# import packages 
import os
import argparse
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
    help='Filename of output.',
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
    help='Test drug category for enrichment.',
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

# split file names and file paths 

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

# set directories
DIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.normpath(os.path.join(DIR, 'DATA'))
OUTDIR = out_path
GENESETDIR = os.path.normpath(os.path.join(DATADIR, 'GENESETS'))
ANNOTDIR = os.path.normpath(os.path.join(DATADIR, 'MAGMA_ANNOT'))

# set filepaths and minimum gene sets size if gene's are named using ENTREZ
if args.id == 'entrez':

    if args.conditional == 'no':
    
        # set gene sets filepaths if setsize is default 
        if args.setsize == 2:
            solo = os.path.normpath(os.path.join(GENESETDIR, 'entrez_genesets.txt'))
            atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_entrez_sets.txt'))
            moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_entrez_sets.txt'))
            ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_entrez_sets.txt'))

        # set file paths for custom minimum gene set size 
        else:
            # create new gene set file for individual drug gene sets 
            df.setsize(GENESETDIR,'/entrez_genesets.txt', args.setsize)
            solo = os.path.normpath(os.path.join(GENESETDIR, 'tmp/entrez_genesets_min%d.txt' % args.setsize))

            # create new gene set file for ATC III code gene sets 
            df.setsize(GENESETDIR,'/atc_entrez_sets.txt', args.setsize)
            atc = os.path.normpath(os.path.join(GENESETDIR, 'tmp/atc_entrez_sets_min%d.txt' % args.setsize))

            # create new gene set file for MOA gene sets 
            df.setsize(GENESETDIR,'/moa_entrez_sets.txt', args.setsize)
            moa = os.path.normpath(os.path.join(GENESETDIR, 'tmp/moa_entrez_sets_min%d.txt' % args.setsize))

            # create new gene set file for clinical indication gene sets 
            df.setsize(GENESETDIR,'/ind_entrez_sets.txt', args.setsize)
            ind = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ind_entrez_sets_min%d.txt' % args.setsize))

    elif args.conditional == 'yes':

        # set gene sets filepaths if setsize is default 
        if args.setsize == 2:
            solo = os.path.normpath(os.path.join(GENESETDIR, 'entrez_cond_sets.txt'))
            atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_cond_sets.txt'))
            moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_cond_sets.txt'))
            ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_cond_sets.txt'))

        # set file paths for custom minimum gene set size 
        else:
            # create new gene set file for individual drug gene sets 
            df.setsize(GENESETDIR,'/entrez_cond_sets.txt', args.setsize)
            solo = os.path.normpath(os.path.join(GENESETDIR, 'tmp/entrez_cond_sets_min%d.txt' % args.setsize))

            # create new gene set file for ATC III code gene sets 
            df.setsize(GENESETDIR,'/atc_cond_sets.txt', args.setsize)
            atc = os.path.normpath(os.path.join(GENESETDIR, 'tmp/atc_cond_sets_min%d.txt' % args.setsize))

            # create new gene set file for MOA gene sets 
            df.setsize(GENESETDIR,'/moa_cond_sets.txt', args.setsize)
            moa = os.path.normpath(os.path.join(GENESETDIR, 'tmp/moa_cond_sets_min%d.txt' % args.setsize))

            # create new gene set file for clinical indication gene sets 
            df.setsize(GENESETDIR,'/ind_cond_sets.txt', args.setsize)
            ind = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ind_cond_sets_min%d.txt' % args.setsize))

# set filepaths and minimum gene sets size if gene's are named using ENSEMBL
elif args.id == 'ensembl':
    
    if args.conditional == 'no':

        # set gene sets filepaths if setsize is default 
        if args.setsize == 2:
            solo = os.path.normpath(os.path.join(GENESETDIR, 'ensembl_genesets.txt'))
            atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_ensembl_sets.txt'))
            moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_ensembl_sets.txt'))
            ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_ensembl_sets.txt'))

        # set file paths for custom minimum gene set size 
        else:
            
            df.setsize(GENESETDIR,'/ensembl_genesets.txt',args.setsize) # individual drug gene sets 
            solo = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ensembl__genesets_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/atc_ensembl_sets.txt',args.setsize) # ATC code gene sets 
            atc = os.path.normpath(os.path.join(GENESETDIR, 'tmp/atc_ensembl_sets_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/moa_ensembl_sets.txt',args.setsize) # MOA gene sets 
            moa = os.path.normpath(os.path.join(GENESETDIR, 'tmp/moa_ensembl_sets_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/ind_ensembl_sets.txt',args.setsize) # clinical indication gene sets 
            ind = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ind_ensembl_sets_min%d.txt' % args.setsize))

    elif args.conditional =='yes':

        # set gene sets filepaths if setsize is default 
        if args.setsize == 2:
            solo = os.path.normpath(os.path.join(GENESETDIR, 'ensembl_cond_sets.txt'))
            atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_ensembl_cond_sets.txt'))
            moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_ensembl_cond_sets.txt'))
            ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_ensembl_cond_sets.txt'))

        # set file paths for custom minimum gene set size 
        else:
            
            df.setsize(GENESETDIR,'/ensembl_cond_sets.txt',args.setsize) # individual drug gene sets 
            solo = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ensembl__cond_sets_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/atc_ensembl_cond_sets.txt',args.setsize) # ATC code gene sets 
            atc = os.path.normpath(os.path.join(GENESETDIR, 'tmp/atc_ensembl_cond_sets_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/moa_ensembl_cond_sets.txt',args.setsize) # MOA gene sets 
            moa = os.path.normpath(os.path.join(GENESETDIR, 'tmp/moa_ensembl_cond_sets_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/ind_ensembl_cond_sets.txt',args.setsize) # clinical indication gene sets 
            ind = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ind_ensembl_cond_sets_min%d.txt' % args.setsize))

# set filepaths and minimum gene sets size if gene's are named using ENSEMBL
elif args.id == 'ensembl92':
    
    if args.conditional == 'no':

        # set gene sets filepaths if setsize is default 
        if args.setsize == 2:
            solo = os.path.normpath(os.path.join(GENESETDIR, 'ensembl_genesets92.txt'))
            atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_ensembl_sets92.txt'))
            moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_ensembl_sets92.txt'))
            ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_ensembl_sets92.txt'))

        # set file paths for custom minimum gene set size 
        else:
            
            df.setsize(GENESETDIR,'/ensembl_genesets92.txt',args.setsize)
            solo = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ensembl_genesets92_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/atc_ensembl_sets92.txt',args.setsize)
            atc = os.path.normpath(os.path.join(GENESETDIR, 'tmp/atcs_ensembl_sets92_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/moa_ensembl_sets92.txt',args.setsize)
            moa = os.path.normpath(os.path.join(GENESETDIR, 'tmp/moa_ensembl_sets92_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/ind_ensembl_sets92.txt',args.setsize)
            ind = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ind_ensembl_sets92_min%d.txt' % args.setsize))

    elif args.conditional == 'yes':

        # set gene sets filepaths if setsize is default 
        if args.setsize == 2:
            solo = os.path.normpath(os.path.join(GENESETDIR, 'ensembl_cond_sets92.txt'))
            atc = os.path.normpath(os.path.join(GENESETDIR, 'atc_ensembl_cond_sets92.txt'))
            moa = os.path.normpath(os.path.join(GENESETDIR, 'moa_ensembl_cond_sets92.txt'))
            ind = os.path.normpath(os.path.join(GENESETDIR, 'ind_ensembl_cond_sets92.txt'))

        # set file paths for custom minimum gene set size 
        else:
            
            df.setsize(GENESETDIR,'/ensembl_cond_sets92.txt',args.setsize)
            solo = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ensembl_cond_sets92_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/atc_ensembl_cond_sets92.txt',args.setsize)
            atc = os.path.normpath(os.path.join(GENESETDIR, 'tmp/atc_ensembl_cond_sets92_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/moa_ensembl_cond_sets92.txt',args.setsize)
            moa = os.path.normpath(os.path.join(GENESETDIR, 'tmp/moa_ensembl_cond_sets92_min%d.txt' % args.setsize))

            df.setsize(GENESETDIR,'/ind_ensembl_cond_sets92.txt',args.setsize)
            ind = os.path.normpath(os.path.join(GENESETDIR, 'tmp/ind_ensembl_cond_sets92_min%d.txt' % args.setsize))


# set MAGMA annotation filepath
if '\\' in args.geneassoc:
    print(not'/' in args.geneassoc or not '\\' in args.geneassoc)

if './' in args.geneassoc or '.\\' in args.geneassoc:
    annot = args.geneassoc.replace('.', os.getcwd(),1)
elif '/' in args.geneassoc or '\\' in args.geneassoc:                        # if no file path is given for the .genes.raw file it assumes it is in the working directory
    annot = args.geneassoc
else:
    annot = os.path.normpath(os.getcwd() + '/' + args.geneassoc)
    
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

# ----------------------------------------------- #
##### ----- PART 4: DRUG GROUP ANALYSIS ----- #####
# ----------------------------------------------- #

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
        df.run_task_silent('Rscript --vanilla %s %s %s %s %s %s' % (os.path.normpath(working+'compute_corrs.R'), annot, os.path.normpath(solo), (gsa), os.path.normpath('/'+out_name), OUTDIR)) 

        # define filepath to set.corrs.rdata and to metadata.rdata file
        corrdata = os.path.normpath(out_path+'/'+'%s_setcorrs.rdata' % out_name)
        metaRdata = os.path.normpath(working+'DATA/metadata.rdata')

        # run either one group or all 

        if args.enrich != 'all':
            
            # compute dependent linear regression
            print('\tRunning dependent linear regression model...')
            df.run_task_silent('Rscript --vanilla %s %s %s %s %s %s %s %s' % (os.path.normpath(working+'compute_lnreg.R'), corrdata, metaRdata, args.enrich.lower(), args.nsize, out_name, (OUTDIR+os.path.sep), args.correct))

        elif args.enrich == 'all':
            
            # define groups to loop through 
            groups = ['atc','moa','ind']

            # loop though types of drug groups and run dependent linear regression model for each type
            for g in groups:

                # compute dependent linear regression
                print('\tRunning dependent linear regression model for %s groups...' % g.upper())
                df.run_task_silent('Rscript --vanilla %s %s %s %s %s %s %s %s' % (os.path.normpath(working+'compute_lnreg.R'), corrdata, metaRdata, g, args.nsize, (out_name+'_'+g.upper()), (OUTDIR+os.path.sep), args.correct))

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


            