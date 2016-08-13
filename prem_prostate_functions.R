
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

	# Add group
	design_df$group <- paste(design_df$cell_line, design_df$days, sep='--')

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

#######################################################
#######################################################
########## 3. Differential expression
#######################################################
#######################################################

#############################################
########## 3.1 Run Voom
#############################################

run_voom <- function(infiles, outfile)
{
	# Load libraries
	library(limma)
	library(edgeR)

	# Load infiles
	design_df <- tread(infiles[1])
	rawcount_mat <- as.matrix(tread(infiles[2], fixRowNames='gene_symbol', subColNames=c('.','-')))

	# Remove DMSO samples
	design_df <- design_df[!grepl('DMSO', design_df$drug),]

	# Define group label
	design_df$group <- paste(design_df$cell_line, design_df$days, sep='--')

	# Replace - with _
	design_df$group <- gsub('-', '_', design_df$group)

	# Select groups to analyze
	selected_groups <- c('LNCaP__10_DAYS', 'LNCaP_RESIDUAL__10_DAYS', 'LNCaP_RESIDUAL__27_DAYS', 'LNCaP_R_CLONES__27_DAYS')

	# Get wells corresponding to the selected groups
	indexes <- which(sapply(design_df$group, function(x) x %in% selected_groups))
	selected_wells <- design_df[indexes, 'well']

	# Get subset of rawcount dataframe
	rawcount_mat <- rawcount_mat[,selected_wells]

	# Make design matrix
	design <- model.matrix(~ 0+names(indexes))

	# Fix column names and replace - with . (- symbols will give errors when running differential expression)
	colnames(design) <- sapply(colnames(design), function(x) gsub('-', '.', strsplit(x, ')')[[1]][2], fixed=TRUE))

	# Create edgeR object
	dge <- DGEList(counts=rawcount_mat)

	# Calculate sample normalization factors
	dge <- calcNormFactors(dge)

	# Normalize with Voom
	v <- voom(dge, design, plot=FALSE)

	# Fit linear model
	fit <- lmFit(v, design)

	# Save results
	save(v, fit, file=outfile)
}

#############################################
########## 3.2 Run differential expression
#############################################

run_differential_expression <- function(infile, outfile)
{
	# Load libraries
	library(limma)

	# Load infile
	load(infile)

	# Get comparison and replace - with _
	comparison <- gsub('-', '_', strsplit(basename(substr(outfile, 1, nchar(outfile)-4)), 'v')[[1]])

	# Prepare contrast string
	contrast_string <- paste(comparison, collapse='-')

	# Make contrast
	contrast <- makeContrasts(contrasts=contrast_string, levels=fit$design)

	# Run differential expression
	fit2 <- contrasts.fit(fit, contrast)
	fit2 <- eBayes(fit2)

	# Get differential expression result table
	de_results_table <- topTable(fit2, n=nrow(fit2$t))

	# Add bonferroni correction
	de_results_table$bonferroni.P.Val <- p.adjust(de_results_table$P.Value, method='bonferroni')

	# Save to outfile
	save(de_results_table, file=outfile)
}

#############################################
########## 3.3 Get logFC table
#############################################

get_logfc_table <- function(infiles, outfile)
{
	# Load infiles
	de_results_list <- lapply(infiles, function(x){ load(x);
											        comparison <- basename(substr(x, 1, nchar(x)-4));
											        de_results_table$comparison <- comparison;
											        de_results_table <- cbind(gene_id=rownames(de_results_table), de_results_table)
											        return(de_results_table)})

	# Make table
	de_results <- do.call('rbind', de_results_list)

	# Get adjusted logFC
	# de_results$logFC_adjusted <- de_results$logFC * (de_results$adj.P.Val < 0.05)*1

	# Get logFC table
	logfc <- reshape2::dcast(gene_id ~ comparison, data=de_results, value.var='logFC')

	# Fix rownames
	rownames(logfc) <- logfc$gene_id
	logfc$gene_id <-  NULL

	# Save results
	save(logfc, file=outfile)
}

#######################################################
#######################################################
########## 4. Coexpression networks
#######################################################
#######################################################

#############################################
########## 4.1 Get coexpression matrices
#############################################

get_coexpression_networks <- function(infiles, outfile)
{
	# Get infiles
	design_df <- tread(infiles[1])
	load(infiles[2])

	# Get condition
	cell_line <- strsplit(basename(outfile), '_coexpression', fixed=TRUE)[[1]][1]

	# Remove cells with timepoints less than 27 days
	design_df <- design_df[design_df$cell_line != 'DU145' & design_df$group != 'LNCaP_RESIDUAL--10_DAYS',]

	# Get wells
	wells <- design_df[design_df$cell_line == cell_line, 'well']

	# Get subset
	vsd_subset <- vsd_df[,wells]

	# Get gene variance
	var <- apply(vsd_subset, 1, var)

	# Remove genes with 0 variance
	vsd_subset <- vsd_subset[names(var)[var > 0],]

	# Get genes to calculate network with
	genes <- c('AR')

	# Calculate correlations
	correlation_list <- lapply(genes,
								function(x)
									lapply(rownames(vsd_subset),
										function(y){
											cor_res <- cor.test(vsd_subset[x,], vsd_subset[y,], method='spearman');
											results <- c(x, y, cor_res$estimate[[1]], cor_res$p.value[[1]])
											return( results ) }))

	# Convert to dataframe
	correlation_network <- as.data.frame(do.call('rbind', do.call('rbind', correlation_list)))

	# Add column names
	colnames(correlation_network) <- c('gene1', 'gene2', 'SCC', 'pvalue')

	# Convert factors to numeric
	correlation_network$SCC <- as.numeric(as.character(correlation_network$SCC))
	correlation_network$pvalue <- as.numeric(as.character(correlation_network$pvalue))

	# Add FDR
	correlation_network$FDR <- p.adjust(correlation_network$pvalue, method='BH')

	# Save result
	save(correlation_network, file=outfile)
}

#######################################################
#######################################################
########## 5. Pathway Enrichment Analysis
#######################################################
#######################################################

#############################################
########## 5.1 Get Pathway enrichment
#############################################

get_pathway_enrichment <- function(infile, outfile)
{
	# Load libraries
	library(citrus)

	# Load infile
	load(infile)

	# Get table type
	table_type <- strsplit(basename(infile), '[_-]+')[[1]][3]

	# Get appropriate table for logFC or msVIPER
	switch(table_type,
		   'msviper' = {gene_score_table <- nes},
		   'logfc' = {gene_score_table <- logfc})

	# Convert rownames
	rownames(gene_score_table) <- s2e(rownames(gene_score_table))

	# Get enrichment
	gsea_nes_table <- apply(gene_score_table, 2, compute_pathway_gsea)

	# Save result
	save(gsea_nes_table, file=outfile)
}

#############################################
########## 5.2 Get Network enrichment
#############################################

get_network_enrichment <- function(infiles, outfile)
{
	# Load libraries
	library(citrus)

	# Get infile labels
	names(infiles) <- gsub('_coexpression.rda', '', basename(infiles))

	# Get coexpression networks
	network_list <- lapply(names(infiles), function(x) {load(infiles[x]);
													    correlation_network$gene2 <- s2e(as.character(correlation_network$gene2));
													    correlation_network$cell_line <- x;
													    return(correlation_network[,c('gene2','SCC','cell_line')])})

	# Concatenate them
	networks <- do.call('rbind', network_list)

	# Cast the table
	correlation_table <- dcast(gene2 ~ cell_line, data=networks, value.var='SCC', fill=0)

	# Fix rownames
	rownames(correlation_table) <- correlation_table$gene2
	correlation_table$gene2 <- NULL
	correlation_table <- as.matrix(correlation_table)

	# Get enrichment
	gsea_nes_table <- apply(correlation_table, 2, compute_pathway_gsea)

	# Save
	save(gsea_nes_table, file=outfile)
}

#############################################
########## 5.3 Make tables
#############################################

get_enrichment_tables <- function(infile, outfile)
{
	# Load infile
	load(infile)

	# Write to table
	twrite(gsea_nes_table, file=outfile, rownames='pathway')
}

#######################################################
#######################################################
########## 6. Modulator Analysis
#######################################################
#######################################################

#############################################
########## 6.1 Get CINDy predictions
#############################################

#############################################
########## 6.2 Calculate Null Model
#############################################

get_resistance_modulator_null <- function(infiles, outfile)
{
	# Load libraries
	library(citrus)

	# Load infiles
	load(infiles[1])
	load(infiles[2])

	# Get prediction list
	prediction <- ifelse(grepl('cindy', infiles[1]), 'cindy', 'preppi')
	switch(prediction,
		   'cindy' = {prediction_list <- prad_cindy_predictions},
		   'preppi' = {prediction_list <- preppi_list})

	# Get NES scores
	nes_scores <- nes[,'LNCaP_R-CLONES--27_DAYSvLNCaP_RESIDUAL--27_DAYS']

	# Convert names to Entrez IDs
	names(nes_scores) <- s2e(names(nes_scores))

	# Get resistance genes
	resistance_genes <- names(nes_scores)[nes_scores > p2z(0.01)]

	# Get null pool
	null_pool <- names(nes_scores)[abs(nes_scores) < p2z(0.5)]

	# Get real counts
	real_counts <- sapply(prediction_list, function(x) length(intersect(x, resistance_genes)))

	# Get N
	n <- as.numeric(strsplit(strsplit(outfile, '-')[[1]][4], '.', fixed=TRUE)[[1]][1])

	# Get null counts
	null_counts <- sapply(paste0('random_', 1:n),
						  function(x) sapply(prediction_list,
						  					 function(y) length(intersect(y, sample(null_pool, length(resistance_genes))))))

	# Make final table
	null_table <- cbind(real=real_counts, null_counts)

	# Save table
	save(null_table, file=outfile)
}

#############################################
########## 6.3 Get Top Modulators
#############################################

get_top_modulators <- function(infile, outfile)
{
	# Load libraries
	library(citrus)

	# Load data
	load(infile)

	# Calculate p-values for each modulator
	pvalues <- apply(null_table, 1, function(x) (1+sum(x[-1] >= x[1])) / (length(x) + 1))

	# Get median number of associated MRs in null
	medians <- apply(null_table[,-1], 1, median)

	# Get result table
	modulator_table <- data.frame(entrez_id = rownames(null_table),
								  gene_symbol = e2s(rownames(null_table)),
								  associated_resistance_genes = null_table[,1],
								  median_associated_random = medians,
								  pvalue = pvalues,
								  FDR = p.adjust(pvalues, method='BH'))

	# Sort
	modulator_table <- modulator_table[order(modulator_table$pvalue, -modulator_table$associated_resistance_genes),]

	# Save table
	twrite(modulator_table, outfile, sort='pvalue')
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

