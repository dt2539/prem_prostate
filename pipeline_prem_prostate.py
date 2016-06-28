
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
########## S1. Data processing
#######################################################
#######################################################

#############################################
########## 1. Make design table
#############################################

@follows(mkdir('s1-expression_tables.dir'))
@files('s0-raw_data.dir/01212015_RNA_ISolation_LIST.xlsx',
	   's1-expression_tables.dir/design_table.txt')

def makeDesignTable(infile, outfile):

	run('make_design_table', infile, outfile, printfiles=False)

##############################
##### 1.1 Plate design plot
##############################

@transform(makeDesignTable,
		  suffix('.txt'),
		  '.png')

def plotPlateDesign(infile, outfile):

	run('plot_plate_design', infile, outfile, printfiles=True)

#############################################
########## 2. Create labeled raw count table
#############################################

@merge((makeDesignTable,
	    's0-raw_data.dir/PE083-A5_R2.cts.txt',
	    's0-raw_data.dir/PE083-A6_R2.cts.txt'),
  		's1-expression_tables.dir/prem_prostate-rawcounts.txt')

def makeRawcountTable(infiles, outfile):

	run('make_rawcount_table', infiles, outfile, printfiles=False, qsub=True)

#############################################
########## 3. Normalize rawcount table
#############################################

@transform(makeRawcountTable,
	       suffix('-rawcounts.txt'),
	       '-normcounts.txt')

def normalizeRawcountTable(infile, outfile):

	run('normalize_rawcount_table', infile, outfile, qsub=True, printfiles=False)

#############################################
########## 4. VST rawcount table
#############################################

@transform(makeRawcountTable,
	       suffix('-rawcounts.txt'),
	       '-vst.txt')

def vstRawcountTable(infile, outfile):

	run('vst_rawcount_table', infile, outfile, qsub=False, printfiles=True)

#############################################
########## 5. Get library sizes
#############################################

@transform(makeRawcountTable,
	       suffix('-rawcounts.txt'),
	       '-library_sizes.txt')

def getLibrarySizes(infile, outfile):

	run('get_library_sizes', infile, outfile, qsub=False, printfiles=False)

##############################
##### 5.1 Histogram plot
##############################

@transform(getLibrarySizes,
		   suffix('.txt'),
		   '.png')

def plotLibrarysizeHistogram(infile, outfile):

	run('plot_librarysize_histogram', infile, outfile, printfiles=False, qsub=False)

##############################
##### 5.2 Plate plot
##############################

@transform(getLibrarySizes,
		   suffix('.txt'),
		   '_plate.png')

def plotLibrarysizePlate(infile, outfile):

	run('plot_librarysize_plate', infile, outfile, printfiles=True)

#######################################################
#######################################################
########## S2. Differential Expression Analysis
#######################################################
#######################################################

#############################################
########## 1. Run Voom
#############################################

@files((makeRawcountTable, makeDesignTable),
	   's2-differential_expression.dir/prem_prostate-voom.rda')

def runVoom(infiles, outfile):

	run('run_voom', infiles, outfile, makedir=True, printfiles=False, qsub=False)

#############################################
########## 2. Run differential expression
#############################################

@transform(runVoom,
		   suffix('voom.rda'),
		   'de_results.rda')

def runDifferentialExpression(infiles, outfile):

	run('run_differential_expression', infiles, outfile, printfiles=False, qsub=False)

#############################################
########## 3. Get resistance DE genes
#############################################

@transform(runDifferentialExpression,
		   suffix('de_results.rda'),
		   'resistance_de_genes.txt')

def getResistanceDeGenes(infile, outfile):

	run('get_resistance_de_genes', infile, outfile, printfiles=False, qsub=False)

##############################
##### 5.1 Venn diagram
##############################

@transform(getResistanceDeGenes,
		   suffix('.txt'),
		   '_venn.png')

def plotDeVenn(infile, outfile):

	run('plot_de_venn', infile, outfile, printfiles=False, qsub=False)

##############################
##### 5.2 Heatmap
##############################

@transform(getResistanceDeGenes,
		   suffix('.txt'),
		   add_inputs(runVoom, makeDesignTable),
		   '_heatmap.png')

def plotDeHeatmap(infile, outfile):

	run('plot_de_heatmap', infile, outfile, printfiles=False, qsub=False)

#######################################################
#######################################################
########## S3. VIPER Analysis
#######################################################
#######################################################

#############################################
########## 1. Merge regulons
#############################################

@merge(glob.glob('/ifs/scratch/c2b2/ac_lab/fmg2117/projects/pancancer/tumAcros005/su2c/*-regulon.rda'),
	   's3-viper.dir/merged_regulon.rda')-

def mergeRegulons(infiles, outfile):

	run('merge_regulons', infiles, outfile, makedir=True, printfiles=False)

#############################################
########## 2. Run msVIPER
#############################################

@files((vstRawcountTable, makeDesignTable, mergeRegulons),
	   's3-viper.dir/prem-prostate_msviper.rda1')

def runMsViper(infiles, outfile):

	run('run_msviper', infiles, outfile, printfiles=True, qsub=False)

#############################################
########## 3. Get resistance MRs
#############################################

@transform(runMsViper,
		   suffix('msviper.rda'),
		   'resistance_mrs.txt')

def getResistanceMrs(infile, outfile):

	run('get_resistance_mrs', infile, outfile, printfiles=False, qsub=False)

##############################
##### 3.1 Venn diagram
##############################

@transform(getResistanceMrs,
		   suffix('.txt'),
		   '_venn.png')

def plotMrVenn(infile, outfile):

	run('plot_mr_venn', infile, outfile, printfiles=False, qsub=False)

#############################################
########## 4. Run ssVIPER
#############################################

@files((vstRawcountTable, makeDesignTable, mergeRegulons),
	   's3-viper.dir/prem_prostate-vipermat.txt')

def runViper(infiles, outfile):

	run('run_viper', infiles, outfile, printfiles=False, qsub=False)

##############################
##### 4.2 Heatmap
##############################

@transform(getResistanceMrs,
		   suffix('.txt'),
		   add_inputs(runViper, makeDesignTable),
		   '_heatmap.png')

def plotMrHeatmap(infile, outfile):

	run('plot_mr_heatmap', infile, outfile, printfiles=False, qsub=False)

#######################################################
#######################################################
########## S4. SU2C Network Analysis
#######################################################
#######################################################

#############################################
########## 1. Network connectivity
#############################################

@files(mergeRegulons,
	   's4-su2c_network.dir/su2c_edge_table.txt')

def getNetworkTables(infile, outfile):

	run('get_network_tables', infile, outfile, makedir=True, printfiles=False)

##############################
##### 1.1 Connectivity histogram
##############################

@transform(getNetworkTables,
		   suffix('_edge_table.txt'),
		   'node_connectivity.png')

def plotNodeConnectivity(infile, outfile):

	run('plot_node_connectivity', infile, outfile, printfiles=False)

##############################
##### 1.2 Edge subset
##############################

@transform(getNetworkTables,
		   suffix('table.txt'),
		   'subset.txt')

def getNetworkSubset(infile, outfile):

	run('get_network_subset', infile, outfile, printfiles=False, qsub=False)





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

