
#################################################################
#################################################################
############### Prostate Cancer Cell Line Analysis - Prem - Functions
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. Custom libraries #####
# Pipeline support
source('/ifs/data/c2b2/ac_lab/dt2539/projects/pipelines/prem_prostate/pipeline/prem_prostate_support.R')

#######################################################
#######################################################
########## S1. Data
#######################################################
#######################################################

#############################################
########## 1. Make design table
#############################################

make_design_table <- function(infile, outfile)
{
	# Load libraries test secondo test
	require(gdata)

	# Read infiles
	platemap_df <- read.xls(infile, sheet='PLATE MAP')
	platemap_df <- as.data.frame(apply(platemap_df, 2, function(x) gsub(' ', '_', x)))

	# Convert days to capitals
	platemap_df$DAYS <- toupper(platemap_df$DAYS)

	# Make design dataframe
	design_df <- data.frame(well = apply(platemap_df[,c('PLATESEQ.PLATE','WELL')], 1, function(x) paste(x, collapse='-')),
							condition = apply(platemap_df[,c('CELL.LINE','TREATMENT','DAYS')], 1, function(x) paste(x, collapse='-')),
							plate = platemap_df[,'PLATESEQ.PLATE'],
							cell_line = platemap_df[,'CELL.LINE'],
							treatment = platemap_df[,'TREATMENT'],
							days = platemap_df[,'DAYS'],
							drug = sapply(as.character(platemap_df[,'TREATMENT']), function(x) strsplit(x, '_')[[1]][1]))

	# Write to outfile
	twrite(design_df, outfile)
}

#############################################
########## 2. Create raw count table
#############################################

make_rawcount_table <- function(infiles, outfile)
{
	# Read infiles
	design_df <- tread(infiles[1])
	rawcount_df_1 <- tread(infiles[2], header=FALSE, fixRowNames='V1')
	rawcount_df_2 <- tread(infiles[3], header=FALSE, fixRowNames='V1')

	# Merge dataframes
	rawcount_df <- cbind(rawcount_df_1, rawcount_df_2)

	# Get column labels
	plate_labels <- design_df$well

	# Remove extra columns
	rawcount_df <- rawcount_df[,1:length(plate_labels)]

	# Remove ERCCs
	ercc_ix <- grepl('ERCC', rownames(rawcount_df))
	rawcount_df <- rawcount_df[!ercc_ix,]

	# Add column labels
	colnames(rawcount_df) <- plate_labels

	# Write to outfile
	twrite(rawcount_df, outfile, rownames='gene_symbol')
}

#############################################
########## 3. VST rawcount table
#############################################

vst_rawcount_table <- function(infile, outfile)
{
	# Load libraries
	library(DESeq2)

	# Read infile
	rawcount_mat <- as.matrix(tread(infile, fixRowNames='gene_symbol', subColNames=c('.','-')))

	# Normalize
	vsd_df <- varianceStabilizingTransformation(rawcount_mat)

	# Write to outfile
	save(vsd_df, file=outfile)
}

#######################################################
#######################################################
########## S2. msVIPER analysis
#######################################################
#######################################################

#############################################
########## 1. Merge SU2C regulons
#############################################

merge_su2c_regulons <- function(infiles, outfile)
{
	# Initialize merged regulon
	regulon <- c()

	# Loop through infiles
	for (infile in infiles)
	{
		# Load regulon
		load(infile)

		# Add
		regulon <- c(regulon, regul)
	}

	# Save to outfile
	save(regulon, file=outfile)
}

#############################################
########## 2. Run msVIPER
#############################################

run_msviper <- function(infiles, outfile)
{	
	# rm(list=ls())
	# source('/ifs/data/c2b2/ac_lab/dt2539/projects/pipelines/prem_prostate/pipeline/prem_prostate_support.R')
	# infiles <-c('f1-data.dir/prem_prostate-vst.rda', 'f1-data.dir/design_table.txt', 'f2-msviper.dir/su2c-merged_regulon.rda')
	# outfile <- 'f2-msviper.dir/prem_prostate-msviper.rda'
	# Load libraries
	library(viper)
	library(citrus)

	# Load data
	load(infiles[1])
	design_df <- tread(infiles[2])
	load(infiles[3])
	
	# Fix rownames	
	rownames(vsd_df) <- s2e(rownames(vsd_df))

	# Remove DMSO samples
	design_df <- design_df[!grepl('DMSO', design_df$drug),]

	# Get groups of interest
	group_labels <- c('LNCaP--10_DAYS', 'LNCaP_RESIDUAL--10_DAYS', 'LNCaP_R-CLONES--27_DAYS', 'LNCaP_RESIDUAL--27_DAYS')

	# Get groups
	groups <- sapply(group_labels, function(x) { cond <- strsplit(x, '--')[[1]];
												 wells <- design_df[(design_df$cell_line == cond[1] & design_df$days == cond[2]), 'well'];
												 return(wells) })

	# Select comparisons
	comparisons <- c('LNCaP--10_DAYSvLNCaP_RESIDUAL--10_DAYS', 'LNCaP_RESIDUAL--10_DAYSvLNCaP_R-CLONES--27_DAYS', 'LNCaP_RESIDUAL--10_DAYSvLNCaP_RESIDUAL--27_DAYS', 'LNCaP_R-CLONES--27_DAYSvLNCaP_RESIDUAL--27_DAYS')
# comparison='LNCaP_R-CLONES--27_DAYSvLNCaP_RESIDUAL--27_DAYS'
	# Define result list
	msviper_results <- list()

	# Run msVIPER
	for (comparison in comparisons)
	{
		# Get pair
		message('\nDoing ', comparison, ' ...')
		pair <- strsplit(comparison, 'v')[[1]]

		# Get samples
		samples1 <- groups[[pair[1]]]
		samples2 <- groups[[pair[2]]]

		# Get matrices
		expmat1 <- vsd_df[,samples1]
		expmat2 <- vsd_df[,samples2]

		# Run t-tests
		sig <- rowTtest(x=expmat1, y=expmat2)
		dnull <- ttestNull(x=expmat1, y=expmat2, per = 1000)

		# Run msVIPER
		good_genes <- rownames(dnull)[complete.cases(dnull)]
		mra <- msviper(sig$statistic[good_genes,], regulon, dnull[good_genes,])
		mrannot <- msviperAnnot(mra, list_eg2symbol)

		# Save to result list
		msviper_results[[comparison]] <- mrannot
	}

	# Save data
	save(msviper_results, file=outfile)
}

#############################################
########## 3. Get NES table
#############################################

get_nes_table <- function(infile, outfile)
{
	# Load infile
	load(infile)

	# Get NES values
	nes_values <- sapply(msviper_results, function(x) x$es$nes)

	# Make melted table list
	nes_table_list <- lapply(names(nes_values), function(x) data.frame(gene_id=names(nes_values[[x]]),
																  NES = nes_values[[x]],
																  comparison = x))

	# Make melted table
	nes_table_merge <- rbind(nes_table_list[[1]], nes_table_list[[2]], nes_table_list[[3]])

	# Cast
	nes_table <- cast(gene_id ~ comparison, data=nes_table_merge, value=c('NES'))

	# Save table
	twrite(nes_table, outfile)
}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################



##################################################
##################################################
########## Run pipeline
##################################################
##################################################
command_args <- commandArgs(trailingOnly = TRUE)
str2R(command_args[1], command_args[2], command_args[3])

