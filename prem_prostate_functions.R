
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
########## S1. Data processing
#######################################################
#######################################################

#############################################
########## 1. Make design table
#############################################

make_design_table <- function(infile, outfile)
{
	# Load libraries test
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

	# Fix cell type for resistant clones
	rclones_ix <- which(design_df$cell_line == 'LNCaP_R-CLONES')
	design_df$cell_type <- as.character(design_df$cell_line)
	design_df[rclones_ix, 'cell_type'] <- paste(design_df[rclones_ix, 'cell_type'], 1:length(rclones_ix), sep='_')

	# Make second condition
	design_df$condition_type <- apply(design_df[,c('cell_type','treatment','days')], 1, function(x) paste(x, collapse='-'))

	# Write to outfile
	twrite(design_df, outfile)
}

##############################
##### 1.1 Plate design plot
##############################

plot_plate_design <- function(infile, outfile)
{
	# Read data
	design_df <- tread(infile)
	design_df$plate <- gsub('_', ' ', design_df$plate)
	head(design_df)

	# Get plate coordinates
	plate_coord_df <- get_table_xy(96, 12)
	plate_coord_df <- rbind(plate_coord_df, plate_coord_df)
	plate_coord_df <- plate_coord_df[1:nrow(design_df),]

	# Concatenate dataframes
	plot_df <- cbind(design_df, plate_coord_df)

	# Get plot ticks and labels
	xticks <- xlabels <- 1:12
	yticks <- 1:8
	ylabels <- toupper(letters)[rev(yticks)]

	# Plot
	png(outfile, height=500, width=1200)

	gp <- ggplot(plot_df, aes(x=x, y=y))
	gp <- gp + geom_point(aes(color=drug), size=15, shape=15, data=plot_df)
	gp <- gp + geom_point(aes(shape=cell_line), size=5, data=plot_df)
	gp <- gp + facet_wrap(~plate, ncol=2)
	gp <- gp + scale_color_manual(values=c('red','turquoise','greenyellow','grey'))
	gp <- gp + scale_x_continuous(breaks=xticks, labels=xlabels)
	gp <- gp + scale_y_continuous(breaks=yticks, labels=ylabels)
	gp <- gp + xlab('') + ylab('') + ggtitle('PLATE-seq experimental design\n')
	print(gp)

	dev.off()

	# # Get colors
	# cell_line_colors <- label_to_color(design_df$cell_line, c('black','slateblue','blue','red'))
	# drug_colors <- label_to_color(design_df$drug, c('turquoise','yellow','green','grey'))
	# plate_colors <- label_to_color(design_df$plate, c('grey60','darkblue'))

}

#############################################
########## 2. Create labeled raw count table
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
########## 3. Normalize rawcount table
#############################################

normalize_rawcount_table <- function(infile, outfile)
{
	# Load libraries
	library(DESeq2)

	# Read infile
	rawcount_df <- tread(infile, fixRowNames='gene_symbol')

	# Calculate size factors
	size_factors <- estimateSizeFactorsForMatrix(rawcount_df)

	# Normalize
	normcount_df <- t(t(rawcount_df/size_factors))

	# Write to outfile
	twrite(normcount_df, outfile, rownames='gene_symbol')
}

#############################################
########## 4. VST rawcount table
#############################################

vst_rawcount_table <- function(infile, outfile)
{
	# Load libraries
	library(DESeq2)
	detach("package:citrus", unload=TRUE)
	detach("package:DESeq", unload=TRUE)

	# Read infile
	rawcount_mat <- as.matrix(tread(infile, fixRowNames='gene_symbol', subColNames=c('.','-')))

	# Normalize
	vsd_df <- varianceStabilizingTransformation(rawcount_mat)

	# Write to outfile
	twrite(vsd_df, outfile, rownames='gene_symbol')
}

#############################################
########## 5. Get library sizes
#############################################

get_library_sizes <- function(infile, outfile)
{
	# Read infile
	rawcount_df <- tread(infile, fixRowNames='gene_symbol', subColNames=c('.','-'))

	# Get column sums
	library_sizes <- colSums(rawcount_df)

	# Make dataframe
	library_size_df <- data.frame(well = names(library_sizes),
								  library_size = library_sizes)

	# Write to outfile
	twrite(library_size_df, outfile)
}

##############################
##### 5.1 Histogram plot
##############################

plot_librarysize_histogram <- function(infile, outfile)
{
	# Read data
	library_size_df <- tread(infile)
	head(library_size_df)

	# Get ticks
	ticks <- seq(0, 2, by=0.5)

	# Plot
	png(outfile, height=800, width=800)

	gp <- ggplot(library_size_df, aes(x=library_size/10^6)) + geom_histogram(bins=100)
	gp <- gp + scale_x_continuous(breaks=ticks, labels=ticks)
	gp <- gp + xlab('Number of reads (x10^6)') + ylab('Count') + ggtitle('PLATE-seq library sizes')
	print(gp)

	dev.off()
}

##############################
##### 5.2 Plate plot
##############################

plot_librarysize_plate <- function(infile, outfile)
{
	# Read data
	library_size_df <- tread(infile)
	library_size_df$plate <- paste('Plate ', sapply(library_size_df$well, function(x) substr(x, 7, 7)))
	head(library_size_df)

	# Get plate coordinates
	plate_coord_df <- get_table_xy(96, 12)
	plate_coord_df <- rbind(plate_coord_df, plate_coord_df)
	plate_coord_df <- plate_coord_df[1:nrow(library_size_df),]

	# Concatenate dataframes
	plot_df <- cbind(library_size_df, plate_coord_df)

	# Get plot ticks and labels
	xticks <- xlabels <- 1:12
	yticks <- 1:8
	ylabels <- toupper(letters)[rev(yticks)]

	# Plot
	png(outfile, height=500, width=1200)

	gp <- ggplot(plot_df, aes(x=x, y=y, color=library_size/10^6)) + geom_point(size=15, shape=15)
	gp <- gp + facet_wrap(~plate, ncol=2)
	gp <- gp + scale_color_gradient2('Library size\n(x10^6)', low='black', mid='yellow', high='red', midpoint=0.8)
	gp <- gp + scale_x_continuous(breaks=xticks, labels=xlabels)
	gp <- gp + scale_y_continuous(breaks=yticks, labels=ylabels)
	gp <- gp + xlab('') + ylab('') + ggtitle('PLATE-seq library sizes\n')
	print(gp)

	dev.off()

}

#######################################################
#######################################################
########## S2. Differential Expression Analysis
#######################################################
#######################################################

#############################################
########## 1. Run celltype voom
#############################################

run_voom <- function(infiles, outfile)
{
	# Load libraries
	library(edgeR)
	library(limma)

	# Read data
	rawcount_df <- tread(infiles[1], fixRowNames='gene_symbol')
	slice(rawcount_df)
	design_df <- tread(infiles[2])
	design_df <- as.data.frame(apply(design_df, 2, function(x) gsub('-', '.', x)))

	# Create design matrix
	design_df$value <- 1
	design <- cast(well ~ cell_line, data=design_df, values='value', fill=0)
	rownames(design) <- design$well
	design$well <- NULL

	# Remove genes with 0 counts
	row_counts <- rowSums(rawcount_df)
	row_ix <- row_counts > 0
	rawcount_df <- rawcount_df[row_ix,]

	# Normalize with VOOM
	dge <- DGEList(counts=rawcount_df)
	dge <- calcNormFactors(dge)
	v <- voom(dge,design,plot=FALSE)

	# Run fit linear model
	fit <- lmFit(v, design)

	# Save data
	save(v, fit, file=outfile)
}

#############################################
########## 2. Run differential expression
#############################################

run_differential_expression <- function(infile, outfile)
{
	# Load libraries
	library(limma)
	library(edgeR)

	# Read data
	load(infile)

	# Get cell pairs
	cell_pairs <- combn(sort(names(v$design)), 2, simplify=FALSE)

	# Initialize result list
	de_list <- list()

	# Loop through cell pairs
	for (cell_pair in cell_pairs)
	{
		# Message
		message('\nDoing ', cell_pair, '...')

		# Get pair label
		pair_label <- gsub('.', '-', paste(cell_pair, collapse='_v_'), fixed=TRUE)

		# Make contrast string
		contrast_string <- paste(cell_pair, collapse='-')
	
		# Get differential expression
		contrast <- makeContrasts(contrasts=contrast_string, levels=fit$design)

		# Fit
		fit2 <- contrasts.fit(fit, contrast)
		fit2 <- eBayes(fit2)

		# Get results
		de_df <- topTable(fit2, n=nrow(fit2$t))

		# Get subset
		de_df <- de_df[,c('logFC','AveExpr','P.Value','adj.P.Val')]

		# Add to list
		de_list[[pair_label]] <- de_df
	}

	# Save to outfile
	save(de_list, file=outfile)
}

#############################################
########## 3. Get resistance DE genes
#############################################

get_resistance_de_genes <- function(infile, outfile)
{
	# Load infile
	load(infile)

	# Get comparison labels
	comparison_labels <- c('LNCaP_v_LNCaP_R-CLONES', 'LNCaP_R-CLONES_v_LNCaP_RESIDUAL')

	# Initialize dataframe
	de_df_melt <- as.data.frame(matrix(nrow=0, ncol=4))

	# Loop through labels
	for (label in comparison_labels)
	{
		# Get DE dataframe
		de_df_temp <- de_list[[label]]
		de_df_temp <- de_df_temp[,c('logFC','adj.P.Val')]
		de_df_temp$gene_id <- rownames(de_df_temp)

		# Melt it
		de_df_temp_melt <- melt(de_df_temp, id=c('gene_id'))

		# Add label
		de_df_temp_melt$id <- paste(de_df_temp_melt$variable, label, sep='_')

		# Append to result
		de_df_melt <- rbind(de_df_melt, de_df_temp_melt)
	}

	# Cast table
	de_df <- cast(gene_id ~ id, data=de_df_melt, value='value')

	# Get combined p-value
	de_df$combined_p <- apply(de_df[,2:3], 1, pstouffer)

	# Get Zscore difference
	de_df$zscore_difference <- apply(de_df[,2:3], 1, function(x) abs( p2z(x[1]) - p2z(x[2]) ))

	# Get combined sign
	de_df$combined_sign <- apply(sign(de_df[,4:5]), 1, prod)

	# Write to outfile
	twrite(de_df, outfile, sort='combined_p')
}

##############################
##### 5.1 Venn diagram
##############################

plot_de_venn <- function(infile, outfile)
{
	# Read data
	de_df <- tread(infile, fixRowNames='gene_id')
	head(de_df)

	# Define list
	de_list <- list()

	# Make lists
	de1 <- rownames(subset(de_df, adj.P.Val_LNCaP_v_LNCaP_R.CLONES < 0.01 & abs(logFC_LNCaP_v_LNCaP_R.CLONES) > 1.5))
	de2 <- rownames(subset(de_df, adj.P.Val_LNCaP_R.CLONES_v_LNCaP_RESIDUAL < 0.01 & abs(logFC_LNCaP_R.CLONES_v_LNCaP_RESIDUAL) > 1.5))
	de_list[['LNCaP_v_LNCaP_R.CLONES']] <- de1
	de_list[['LNCaP_R-CLONES_v_LNCaP_RESIDUAL']] <- de2

	# Plot venn
	png(outfile, height=800, width=800)

	venn(de_list)	

	dev.off()

}

##############################
##### 5.2 Heatmap
##############################

plot_de_heatmap <- function(infiles, outfile)
{
	# Read data
	de_df <- tread(infiles[1], fixRowNames='gene_id')
	head(de_df)
	load(infiles[2])
	design_df <- tread(infiles[3], fixRowNames='well')
	head(design_df)

	# Get normcount dataframe
	normcount_df <- v$E
	colnames(normcount_df) <- gsub('.', '-', colnames(normcount_df), fixed=TRUE)

	# Get data subset
	wells <- rownames(subset(design_df, cell_line != 'DU145'))
	design_df <- design_df[wells,]
	normcount_df <- normcount_df[,wells]

	# Get DE genes
	de1 <- rownames(subset(de_df, adj.P.Val_LNCaP_v_LNCaP_R.CLONES < 0.01 & abs(logFC_LNCaP_v_LNCaP_R.CLONES) > 1))
	de2 <- rownames(subset(de_df, adj.P.Val_LNCaP_R.CLONES_v_LNCaP_RESIDUAL < 0.01 & abs(logFC_LNCaP_R.CLONES_v_LNCaP_RESIDUAL) > 1))
	# de2 <- rownames(subset(de_df, adj.P.Val_LNCaP_v_LNCaP_RESIDUAL < 0.01 & abs(logFC_LNCaP_v_LNCaP_RESIDUAL) > 1))
	de <- intersect(de1, de2)

	# Get plot dataframe
	plot_df <- normcount_df[de,]

	# Prepare colors
	cell_line_colors <- label_to_color(design_df$cell_line)
	drug_colors <- label_to_color(design_df$drug, topo.colors(4))

	# Plot
	png(outfile, height=1000, width=1000)

	heatmap.3(plot_df,
			  ColSideColors = data.frame('Cell type'=cell_line_colors$colors, 'Treatment'=drug_colors$colors),
			  keysize = 1,
			  scale = 'row',
			  KeyValueName = 'Expression level (Voom)',
			  cexRow = 1,
			  margins = c(5, 7),
			  distfun = function(x) as.dist(1-cor(t(x), method='spearman')))

	legend(0.2, 1, names(cell_line_colors$legend), fill=as.character(cell_line_colors$legend), title='Cell type', bg='white')
	legend(0.4, 1, names(drug_colors$legend), fill=as.character(drug_colors$legend), title='Treatment', bg='white')

	dev.off()

}

#######################################################
#######################################################
########## S3. VIPER Analysis
#######################################################
#######################################################

#############################################
########## 1. Merge regulons
#############################################

merge_regulons <- function(infiles, outfile)
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

	# Read data
	normcount_df <- as.matrix(tread(infiles[1], fixRowNames='gene_symbol', subColNames=c('.','-')))
	design_df <- tread(infiles[2])
	load(infiles[3])
	slice(normcount_df)
	head(design_df)

	# Convert names
	rownames(normcount_df) <- s2e(rownames(normcount_df))

	# Get cell2plate list
	cell2plate <- sapply(unique(design_df$cell_line), function(x) design_df[design_df$cell_line == x, 'well'])

	# Get cell pairs
	cell_pairs <- combn(unique(design_df$cell_line), 2, simplify=FALSE)

	# Define results
	msviper_results <- list()

	# Load eg2symbol list
	e2s('10001')

	# Loop through cell pairs
	for (cell_pair in cell_pairs)
	{
		# Message
		message('\nDoing ', cell_pair, '...')

		# Get pair label
		pair_label <- paste(cell_pair, collapse='_v_')

		# Get wells
		wells1 <- cell2plate[[cell_pair[1]]]
		wells2 <- cell2plate[[cell_pair[2]]]

		# Run msVIPER
		sig <- rowTtest(x=normcount_df[,wells1], y=normcount_df[,wells2])
		dnull <- ttestNull(x=normcount_df[,wells1], y=normcount_df[,wells2], per = 1000)
		good_genes <- rownames(dnull)[complete.cases(dnull)]
		mra <- msviper(sig$statistic[good_genes,], regulon, dnull[good_genes,])
		mrannot <- msviperAnnot(mra, list_eg2symbol)

		# mrs <- mra
		# rownames(mrs$signature) <- e2s(rownames(mrs$signature))
		# pdf()
		# plot(mrs)
		# dev.off()

		# Save to list
		msviper_results[[pair_label]] <- mrannot
	}

	# Save data
	save(msviper_results, file=outfile)
}

#############################################
########## 3. Get resistance MRs
#############################################

get_resistance_mrs <- function(infile, outfile)
{
	# Load data
	load(infile)

	# Get comparison names
	comparison_labels <- c('LNCaP_v_LNCaP_R-CLONES', 'LNCaP_RESIDUAL_v_LNCaP_R-CLONES')

	# Initialize list
	nes_df_melt <- as.data.frame(matrix(nrow=0, ncol=3))

	# Loop through labels
	for (label in comparison_labels)
	{
		# Get NES vector
		nes_vector <- msviper_results[[label]]$es$nes
		nes_adjusted <- nes_vector#p2z(p.adjust(z2p(nes_vector), method='BH'))

		# Create dataframe
		nes_df_temp <- data.frame(gene_id = names(nes_vector),
								  nes = nes_adjusted,
								  comparison = label)

		# Append
		nes_df_melt <- rbind(nes_df_melt, nes_df_temp)
	}

	# Cast dataframe
	nes_df <- cast(gene_id ~ comparison, data=nes_df_melt, value='nes')

	# Get combined Z-scores
	nes_df$combined_z <- apply(abs(nes_df[,comparison_labels]), 1, stouffer)

	# Get sign
	nes_df$combined_sign <- apply(sign(nes_df[,comparison_labels]), 1, prod)

	# Get intersection
	p <- 0.01
	print(head(nes_df))
	nes_df$both <- abs(nes_df[,2]) > p2z(p) & abs(nes_df[,3]) > p2z(p) 

	# Save
	twrite(nes_df, outfile, sort='combined_z', decreasing=TRUE)
}

##############################
##### 3.1 Venn diagram
##############################

plot_mr_venn <- function(infile, outfile)
{
	# Read data
	nes_df <- tread(infile, fixRowNames='gene_id')
	head(nes_df)

	# Define list
	mr_list <- list()

	# Make lists
	for (condition in colnames(nes_df)[1:2])
	{
		# Add MRs
		mrs <- rownames(nes_df)[abs(nes_df[,condition]) > p2z(0.01)]
		mr_list[[condition]] <- mrs[!is.na(mrs)]
	}

	# Plot venn
	png(outfile, height=800, width=800)

	venn(mr_list)	

	dev.off()
}

#############################################
########## 4. Run ssVIPER
#############################################

run_viper <- function(infiles, outfile)
{
	# Load libraries
	library(viper)

	# Read data
	normcount_df <- as.matrix(tread(infiles[1], fixRowNames='gene_symbol', subColNames=c('.','-')))
	design_df <- tread(infiles[2])
	load(infiles[3])
	slice(normcount_df)
	head(design_df)

	# Convert names
	rownames(normcount_df) <- s2e(rownames(normcount_df))

	# Run VIPER
	viper_results <- viper(normcount_df, regulon)

	# Save results
	twrite(viper_results, outfile, rownames='gene_id')
}

##############################
##### 4.2 Heatmap
##############################

plot_mr_heatmap <- function(infiles, outfile)
{
	# Read data
	mr_df <- tread(infiles[1], fixRowNames='gene_id')
	head(mr_df)
	viper_df <- tread(infiles[2], fixRowNames='gene_id', subColNames=c('.','-'))
	slice(viper_df)
	design_df <- tread(infiles[3], fixRowNames='well')
	head(design_df)

	# Fix rownames
	rownames(viper_df) <- e2s(rownames(viper_df))

	# Get data subset
	wells <- rownames(subset(design_df, cell_line != 'DU145'))
	design_df <- design_df[wells,]
	viper_df <- viper_df[,wells]

	# Get DE genes
	mr1 <- rownames(subset(mr_df, abs(LNCaP_v_LNCaP_R.CLONES) > p2z(0.01)))
	mr2 <- rownames(subset(mr_df, abs(LNCaP_RESIDUAL_v_LNCaP_R.CLONES) > p2z(0.01)))
	mr <- intersect(mr1, mr2)

	# Get plot dataframe
	plot_df <- viper_df[mr,]

	# Prepare colors
	cell_line_colors <- label_to_color(design_df$cell_line)
	drug_colors <- label_to_color(design_df$drug, topo.colors(5))

	# Plot
	png(outfile, height=1300, width=1300, pointsize=17)

	heatmap.3(plot_df,
			  ColSideColors = data.frame('Cell type'=cell_line_colors$colors, 'Treatment'=drug_colors$colors),
			  keysize = 1,
			  distfun = function(x) as.dist(1-cor(t(x), method='pearson')),
			  cexRow=1,
			  KeyValueName = 'NES')

	legend(0.2, 1, names(cell_line_colors$legend), fill=as.character(cell_line_colors$legend), title='Cell type', bg='white')
	legend(0.4, 1, names(drug_colors$legend), fill=as.character(drug_colors$legend), title='Treatment', bg='white')

	dev.off()

}

#######################################################
#######################################################
########## S4. SU2C Network Analysis
#######################################################
#######################################################

#############################################
########## 1. Network plot
#############################################

get_network_tables <- function(infile, outfile)
{
	# Load infile
	load(infile)

	# Define edge table
	edge_df <- as.data.frame(matrix(nrow=0, ncol=4))

	# Fill table
	for (tf in names(regulon))
	{
		# Make dataframe
		tf_edge_df <- data.frame(TF = tf,
									   target = names(regulon[[tf]]$tfmode),
									   sign = sign(regulon[[tf]]$tfmode),
									   likelihood = regulon[[tf]]$likelihood)

		# Append data
		edge_df <- rbind(edge_df, tf_edge_df)
	}
	
	# Define node table
	tf_df <- data.frame(id = names(regulon), label='TF')
	target_df <- data.frame(id = unique(edge_df$target), label='target')
	node_df <- rbind(tf_df, target_df)
	node_df$symbol <- e2s(node_df$id)

	# Save files
	node_outfile <- gsub('edge', 'node', outfile)
	twrite(edge_df, outfile)
	twrite(node_df, node_outfile)
}

##############################
##### 1.1 Connectivity histogram
##############################

plot_node_connectivity <- function(infile, outfile)
{
	# Read table
	edge_df <- tread(infile)
	head(edge_df)

	# Get TF connectivity
	connectivity_df <- aggregate(target ~ TF, data=edge_df, FUN=length)
	connectivity_df$TF <- e2s(connectivity_df$TF)
	colnames(connectivity_df)[2] <- 'nr_targets'
	head(connectivity_df)

	# Plot
	png(outfile, height=500, width=500)

	gp <- ggplot(connectivity_df, aes(x=nr_targets)) + geom_histogram(bins=200)
	gp <- gp + geom_text(aes(label=TF, y=10, angle=60, hjust=0), data=subset(connectivity_df, nr_targets > 600 & TF != 'EYA4'))
	gp <- gp + scale_x_continuous(lim=c(0, 1000))
	gp <- gp + xlab('Number of targets') + ylab('Count') + ggtitle('Node Connectivity Distribution')
	gp

	dev.off()

}

##############################
##### 1.2 Edge subset
##############################

get_network_subset <- function(infile, outfile)
{
	# Read table
	edge_df <- tread(infile)
	head(edge_df)

	# Get genes
	# gene_symbols <- c('SPDEF', 'NKX3-1', 'WDR77', 'STMN3')
	# gene_symbols <- c('KLF15', 'ZNF613', 'SP4', 'ZNF649', 'ZNF615', 'KDM2B')
	gene_symbols <- c('STAT1', 'ERBB2', 'PTPRJ', 'EIF2AK3', 'IFI6')
	gene_ids <- s2e(gene_symbols)
	print(gene_ids)

	# Get subset
	edge_df_subset <- subset(edge_df, TF %in% gene_ids)

	# Convert to gene symbol
	edge_df_subset$TF <- e2s(edge_df_subset$TF)
	edge_df_subset$target <- e2s(edge_df_subset$target)

	# Write to outfile
	twrite(edge_df_subset, outfile)
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

