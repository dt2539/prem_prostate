
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
msigdbGenesets = '/ifs/data/c2b2/ac_lab/dt2539/projects/project_data/f1-genelists.dir/msigdb-genesets.rda'
cindyPredictions = '/ifs/data/c2b2/ac_lab/dt2539/projects/project_data/f2-interaction_lists.dir/prad_cindy-list.rda'
preppiPredictions = '/ifs/data/c2b2/ac_lab/dt2539/projects/project_data/f2-interaction_lists.dir/preppi-list.rda'

##### 2. Other variables #####
# Set comparisons
comparisons = ['LNCaP--10_DAYSvLNCaP_RESIDUAL--10_DAYS',
			   'LNCaP_RESIDUAL--10_DAYSvLNCaP_R-CLONES--27_DAYS',
			   'LNCaP_RESIDUAL--10_DAYSvLNCaP_RESIDUAL--27_DAYS',
			   'LNCaP_R-CLONES--27_DAYSvLNCaP_RESIDUAL--27_DAYS',
			   'LNCaP--10_DAYSvLNCaP_RESIDUAL--27_DAYS',
			   'LNCaP--10_DAYSvLNCaP_R-CLONES--27_DAYS']

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

def msviperJobs():
	# Set infiles
	infiles = ['f1-data.dir/prem_prostate-vst.rda',
			   'f1-data.dir/design_table.txt',
			   'f2-msviper.dir/su2c-merged_regulon.rda']
	# Loop through comparisons
	for comparison in comparisons:
		outfile = 'f2-msviper.dir/msviper_runs/%(comparison)s.rda' % locals()
		yield [infiles, outfile]

@follows(mkdir('f2-msviper.dir/msviper_runs/'))

@files(msviperJobs)

def runMsViper(infiles, outfile):

	run('run_msviper', infiles, outfile, mem='32G')

#############################################
########## 3. Get NES table
#############################################

### Purpose: Get msVIPER NES results table.

@merge(runMsViper,
	   'f2-msviper.dir/prem_prostate-msviper_table.rda')

def getNesTable(infiles, outfile):

	run('get_nes_table', infiles, outfile)

#############################################
########## 4. Get significant gene list
#############################################

@transform(getNesTable,
		   suffix('table.rda'),
		   'significant_genes.rda')

def getSignificantGenes(infile, outfile):

	run('get_significant_genes', infile, outfile)

#############################################
########## 5. Get significant gene table
#############################################

@transform(getSignificantGenes,
		   suffix('.rda'),
		   '.txt')

def getSignificantGeneTable(infile, outfile):

	run('get_significant_gene_table', infile, outfile)

#######################################################
#######################################################
########## 3. Differential expression
#######################################################
#######################################################

#############################################
########## 3.1 Run edgeR normalization
#############################################

@follows(mkdir('f3-differential_expression.dir/de_runs'))

@files((makeDesignTable, makeRawcountTable),
	   'f3-differential_expression.dir/prem_prostate-voom.rda')

def runVoom(infiles, outfile):

	run('run_voom', infiles, outfile)

#############################################
########## 3.2 Run differential expression
#############################################

def deJobs():
	# Set infiles
	infile = 'f3-differential_expression.dir/prem_prostate-voom.rda'
	# Loop through comparisons
	for comparison in comparisons:
		outfile = 'f3-differential_expression.dir/de_runs/%(comparison)s.rda' % locals()
		yield [infile, outfile]

@follows(mkdir('f3-differential_expression.dir/de_runs'))

@files(deJobs)

def runDifferentialExpression(infile, outfile):

	run('run_differential_expression', infile, outfile)

#############################################
########## 3.3 Get logFC table
#############################################

@files(runDifferentialExpression,
	   'f3-differential_expression.dir/prem_prostate-logfc_table.rda')

def getLogfcTable(infiles, outfile):

	run('get_logfc_table', infiles, outfile)

#######################################################
#######################################################
########## 4. Coexpression networks
#######################################################
#######################################################

#############################################
########## 4.1 Get coexpression matrices
#############################################

def coexpressionJobs():
	# Get infile
	infiles = ['f1-data.dir/design_table.txt', 'f1-data.dir/prem_prostate-vst.rda']
	for cell_line in ['LNCaP', 'LNCaP_RESIDUAL', 'LNCaP_R-CLONES']:
		outfile = 'f4-coexpression.dir/%(cell_line)s_coexpression.rda' % locals()
		yield [infiles, outfile]

@follows(mkdir('f4-coexpression.dir'))

@files(coexpressionJobs)

def getCoexpressionNetworks(infiles, outfile):

	run('get_coexpression_networks', infiles, outfile)

#############################################
########## 4.2 Compare networks
#############################################

@merge(getCoexpressionNetworks,
	   'f4-coexpression.dir/prem_prostate-network_comparison.txt')

def compareCoexpressionNetworks(infiles, outfile):

	run('compare_coexpression_networks', infiles, outfile)

#######################################################
#######################################################
########## 5. Pathway Enrichment Analysis
#######################################################
#######################################################

#############################################
########## 5.1 Get Pathway enrichment
#############################################

@follows(mkdir('f5-pathway_enrichment.dir'))

@transform((getNesTable, getLogfcTable),
		   regex(r'.*/(.*).rda'),
		   r'f5-pathway_enrichment.dir/\1_enrichment.rda')

def getPathwayEnrichment(infile, outfile):

	run('get_pathway_enrichment', infile, outfile, mem='16G')

#############################################
########## 5.2 Get Network enrichment
#############################################

@merge(getCoexpressionNetworks,
	   'f5-pathway_enrichment.dir/prem_prostate-AR_network_enrichment.rda')

def getNetworkEnrichment(infiles, outfile):

	run('get_network_enrichment', infiles, outfile, mem='16G')

#############################################
########## 5.3 Make tables
#############################################

@transform('f5-pathway_enrichment.dir/prem_prostate-logfc_table_enrichment.rda',#(getPathwayEnrichment, getNetworkEnrichment),
		   suffix('.rda'),
		   '.txt')

def getEnrichmentTables(infile, outfile):

	run('get_enrichment_tables', infile, outfile)

#######################################################
#######################################################
########## 6. Modulator Analysis
#######################################################
#######################################################

#############################################
########## 6.1 Calculate Null Model
#############################################

@transform((preppiPredictions, cindyPredictions),
		   regex(r'.*/(.*)-list.rda'),
		   add_inputs(getNesTable),
		   r'f6-modulator_analysis.dir/prem_prostate-resistance_genes_\1_null-10000.rda')

def getResistanceModulatorNull(infiles, outfile):

	run('get_resistance_modulator_null', infiles, outfile)

#############################################
########## 6.2 Get Top Modulators
#############################################

@transform(getResistanceModulatorNull,
		   suffix('.rda'),
		   '_predictions.txt')

def getTopModulators(infile, outfile):

	run('get_top_modulators', infile, outfile)


#######################################################
#######################################################
########## . 
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

