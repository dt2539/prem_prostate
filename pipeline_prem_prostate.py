
#################################################################
#################################################################
############### Prostate Cancer Cell Line Analysis - Prem
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import glob, sys, os

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/ifs/data/c2b2/ac_lab/dt2539/projects/scripts')
from py2str import *

# Pipeline support
sys.path.append('/ifs/data/c2b2/ac_lab/dt2539/projects/support')
import support as S

#############################################
########## 2. R connection
#############################################
##### 1. R Connection setup #####
pipeline_functions_source = '/ifs/data/c2b2/ac_lab/dt2539/projects/pipelines/prem_prostate/pipeline/prem_prostate_functions.R'
working_dir = '/ifs/data/c2b2/ac_lab/dt2539/projects/pipelines/prem_prostate'
def run(function, infiles, outfile, makedir=False, qsub=True, mem='8000M', nodes=1, printfiles=False):
	runR(function, infiles, outfile, makedir, pipeline_functions_source, working_dir, qsub, mem, nodes, printfiles)

#############################################
########## 3. General setup
#############################################
##### 1. Default variables #####


#######################################################
#######################################################
########## S1. Data
#######################################################
#######################################################

#############################################
########## 1. Make design table
#############################################

### Purpose: Create a design table for the PCA PLATE-Seq experiment.

@follows(mkdir('f1-data.dir'))
@files('f0-raw_data.dir/01212015_RNA_ISolation_LIST.xlsx',
	   'f1-data.dir/design_table.txt')

def makeDesignTable(infile, outfile):

	run('make_design_table', infile, outfile)

#############################################
########## 2. Create raw count table
#############################################

### Purpose: Create merged rawcount table.

@merge((makeDesignTable,
	    'f0-raw_data.dir/PE083-A5_R2.cts.txt',
	    'f0-raw_data.dir/PE083-A6_R2.cts.txt'),
  		'f1-data.dir/prem_prostate-rawcounts.txt')

def makeRawcountTable(infiles, outfile):

	run('make_rawcount_table', infiles, outfile)

#############################################
########## 3. VST rawcount table
#############################################

### Purpose: Run VST transformation on the rawcount table.

@transform(makeRawcountTable,
	       suffix('rawcounts.txt'),
	       'vst.rda')

def vstRawcountTable(infile, outfile):

	run('vst_rawcount_table', infile, outfile)

#######################################################
#######################################################
########## S2. msVIPER analysis
#######################################################
#######################################################

#############################################
########## 1. Merge SU2C regulons
#############################################

### Purpose: Merge TF, coTF and signaling regulons from SU2C.

@follows(mkdir('f2-msviper.dir'))

@merge(glob.glob('/ifs/scratch/c2b2/ac_lab/fmg2117/projects/pancancer/tumAcros005/su2c/*-regulon.rda'),
	   'f2-msviper.dir/su2c-merged_regulon.rda')

def mergeSu2cRegulons(infiles, outfile):

	run('merge_su2c_regulons', infiles, outfile)

#############################################
########## 2. Run msVIPER
#############################################

### Purpose: Run msVIPER.

@files((vstRawcountTable, makeDesignTable, mergeSu2cRegulons),
	   'f2-msviper.dir/prem_prostate-msviper.rda')

def runMsViper(infiles, outfile):

	run('run_msviper', infiles, outfile)

#############################################
########## 3. Get NES table
#############################################

### Purpose: Get msVIPER NES results table.

@transform('f2-msviper.dir/prem_prostate-msviper_100.rda',
		   suffix('.rda'),
		   '.txt')

def getNesTable(infile, outfile):

	run('get_nes_table', infile, outfile)


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

##############################
##### .
##############################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')

