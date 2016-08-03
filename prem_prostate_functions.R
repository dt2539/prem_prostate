
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
	group_labels <- strsplit(basename(substr(outfile, 1, nchar(outfile)-4)), 'v')[[1]]

	# Get groups
	groups <- lapply(group_labels, function(x) { cond <- strsplit(x, '--')[[1]];
												 wells <- design_df[(design_df$cell_line == cond[1] & design_df$days == cond[2]), 'well'];
												 return(wells) })

	# Get matrices
	expmats <- lapply(groups, function(x) vsd_df[,x])

	# Run t-tests
	sig <- rowTtest(x=expmats[[1]], y=expmats[[2]])
	dnull <- ttestNull(x=expmats[[1]], y=expmats[[2]], per = 10000)

	# Run msVIPER
	good_genes <- rownames(dnull)[complete.cases(dnull)]
	mra <- msviper(sig$statistic[good_genes,], regulon, dnull[good_genes,])
	msviper_results <- msviperAnnot(mra, list_eg2symbol)

	# Save data
	save(msviper_results, file=outfile)
}

#############################################
########## 3. Get NES table
#############################################

get_nes_table <- function(infiles, outfile)
{
	# Load libraries
	library(citrus)

	# Get NES data
	nes_data <- lapply(infiles, function(x) { load(x);
											  label <- strsplit(basename(x), '.', fixed=TRUE)[[1]][1]
											  nes <- msviper_results$es$nes
											  nes_table <- data.frame(gene_id = as.character(names(nes)),
											  						  nes = nes,
											  						  comparison = as.character(label))
											  return(nes_table)})

	# Get NES table
	nes_table <- do.call('rbind', nes_data)

	# Cast it
	nes <- reshape2::dcast(gene_id ~ comparison, data=nes_table, value.var='nes')

	# Fix row names
	rownames(nes) <- nes$gene_id
	nes$gene_id <- NULL

	# Remove NA
	nes_filter <- as.matrix(nes[complete.cases(nes),])

	# Adjust p-values
	nes <- apply(nes_filter, 2, function(x) p2z(p.adjust(z2p(x), method='BH'))*sign(x))

	# Fix rownames
	rownames(nes) <- rownames(nes_filter)

	# Save table
	save(nes, file=outfile)
}

#############################################
########## 4. Get significant gene list
#############################################

get_significant_genes <- function(infile, outfile)
{
	# Load infiles
	load(infile)

	# Fix sign
	nes <- -nes

	# Set p threshold
	p_threshold <- 0.05

	# Get VIPER genelists
	significant_genes <- list()

	# Loop through comparisons
	for (comparison in colnames(nes))
	{
		for (sign in c('up', 'down'))
		{
			# Set label
			label <- paste(comparison, sign, sep='__')

			# Add results
			significant_genes[[label]] <- switch(sign,
											   'up' = rownames(nes)[nes[,comparison] > citrus::p2z(p_threshold)],
											   'down' = rownames(nes)[nes[,comparison] < -citrus::p2z(p_threshold)])
		}
	}

	# Save results
	save(significant_genes, file=outfile)
}

#############################################
########## 5. Get significant gene table
#############################################

get_significant_gene_table <- function(infile, outfile)
{
	# Load infile
	load(infile)

	# Get max length of sets
	max_length <- max(sapply(significant_genes, length))

	# Get updated list
	significant_gene_table <- sapply(names(significant_genes), function(x) {genes <- significant_genes[[x]];
																		    genes <- c(genes, rep('', max_length-length(genes)));
																		    names(genes) <- rep(x, length(genes));
																		    return(genes) })

	# Save table
	twrite(significant_gene_table, outfile)
}

#############################################
########## 5. Get pathway enrichment
#############################################

get_pathway_enrichment <- function(infiles, outfile)
{
	# Load libraries
	library(citrus)

	# Load infiles
	load(infiles[1])
	load(infiles[2])
	load(infiles[3])

	# Get background
	background <- s2e(rownames(nes))

	# Intersect with background
	msigdb_genesets_filtered <- sapply(msigdb_genesets, function(x) intersect(x, background))

	# Filter by size
	geneset_size <- sapply(msigdb_genesets_filtered, length)
	msigdb_genesets_filtered <- msigdb_genesets_filtered[names(geneset_size)[geneset_size > 30]]

	# Convert gene IDs
	significant_genes <- sapply(significant_genes, function(x) s2e(x))

	# Define fisher test result list
	fisher_results <- list()

	get_pathway_enrichment <- function(geneset, pathway_genes, background)
	{
	    # Get binary table
	    binary_table <- data.frame(in_geneset = sapply(background, function(x) x %in% geneset),
	                               in_pathway = sapply(background, function(x) x %in% pathway_genes))

	    # Get contingency table
	    contingency_table <- table(binary_table)

	    # Run Fisher's test
	    fisher_test_results <- tryCatch(fisher.test(contingency_table), error=function(x) NA)

	    # Return result
	    return(fisher_test_results)
	}



	# Loop through genesets
	for (geneset_name in names(msigdb_genesets_filtered))
	{
		for (comparison_name in names(significant_genes))
		{
			# Get label
			label <- paste(geneset_name, comparison_name, sep='___')

			# Run fisher test
			fisher_results[[label]] <- get_pathway_enrichment(geneset = significant_genes[[comparison_name]],
													 		  pathway_genes = msigdb_genesets_filtered[[geneset_name]],
															  background = background)
		}
	}

	# Get pvalues
	pvalues <- data.frame(p=sapply(fisher_results, function(x) tryCatch(x$p.value, error=function(x) 1)))
	pvalues$hallmark <- sapply(rownames(pvalues), function(x) strsplit(x, '___')[[1]][1])
	pvalues$comparison <- sapply(rownames(pvalues), function(x) strsplit(x, '___')[[1]][2])
	pvalues$FDR <- p.adjust(pvalues$p, method='BH')
	pvalues$z <- p2z(pvalues$p)

	# Cast table
	zscore_table <- reshape2::dcast(comparison ~ hallmark, data=pvalues, value.var='z', fill=1)
	rownames(zscore_table) <- zscore_table$comparison
	zscore_table$comparison <- NULL

	# Save results
	save(zscore_table, fisher_results, file=outfile)
}

#############################################
########## 6. Get resistance genes
#############################################

get_resistance_genes <- function(infiles, outfile)
{
	# Load infile
	load(infiles[1])
	load(infiles[2])

	# Get resistance genes
	resistance_genes <- setdiff(significant_genes[['LNCaP_R-CLONES--27_DAYSvLNCaP_RESIDUAL--27_DAYS__down']],
								significant_genes[['LNCaP--10_DAYSvLNCaP_RESIDUAL--10_DAYS__down']])

	# Get NES
	resistance_gene_nes <- sort(nes[resistance_genes, 'LNCaP_R-CLONES--27_DAYSvLNCaP_RESIDUAL--27_DAYS'])



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

