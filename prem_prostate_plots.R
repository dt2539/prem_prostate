
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
########## S2. msVIPER analysis
#######################################################
#######################################################

#############################################
########## 2.1 Cell cycle scatter
#############################################

### Input and output
infile <- 'f2-msviper.dir/prem_prostate-msviper_100.txt'
outfile <- 'plots/prem_prostate-msviper_scatter.png'

### Code
nes_table <- tread(infile)

# Get colors
cc_table <- tread('plots/cc_genes.txt')
cc_genes <- e2s(unique(cc_table[,1]))
cc_binary <- sapply(nes_table$gene_id, function(x) x %in% cc_genes)
cc_colors <- sapply(cc_binary, function(x) ifelse(x, 'red', 'black'))

# Open device
png(outfile, height=500, width=500)

plot(-nes_table$LNCaPvLNCaP_RESIDUAL, -nes_table$LNCaPvLNCaP_R.CLONES, col=cc_colors,
	 xlab='NES\nLNCaP Untreated to Residual', ylab='NES\nLNCaP Residual to Resistant',
	 main='msVIPER analysis\nTransition between LNCaP Untreated, Residual, Resistant',
	 xlim = c(-4,4), ylim=c(-4,4))

p <- 0.05
abline(h=p2z(p), lty=2)
abline(h=-p2z(p), lty=2)
abline(v=p2z(p), lty=2)
abline(v=-p2z(p), lty=2)
abline(a=0, b=1)

dev.off()

#############################################
########## 2.2 Scatter plot
#############################################

# Make directory
dir.create('plots/f2-msviper')

### Input and output
infile <- 'f2-msviper.dir/prem_prostate-msviper_table.rda'
outfile <- 'plots/f2-msviper/prem_prostate-msviper_scatter.png'

### Code
# Load data
load(infile)
nes <- nes[complete.cases(nes),]

# Make color map
nes.cmap <- squash::makecmap(nes[,'LNCaP--10_DAYSvLNCaP_RESIDUAL--10_DAYS'], n=30, colFn=colorRampPalette(c('navy', 'white', 'red3')))

# Find interesting genes
genes <- c()
genes <- c(genes, rownames(nes[nes[,'LNCaP--10_DAYSvLNCaP_RESIDUAL--27_DAYS'] < -2 & nes[,'LNCaP--10_DAYSvLNCaP_R-CLONES--27_DAYS'] > 2,]))
genes <- c(genes, rownames(nes[nes[,'LNCaP--10_DAYSvLNCaP_RESIDUAL--27_DAYS'] > 2 & nes[,'LNCaP--10_DAYSvLNCaP_R-CLONES--27_DAYS'] < -2,]))
genes <- c(genes, c('RASD2', 'IRF6', 'NFKB2'))

# Plot scatter
png(outfile, height=1400, width=1600, pointsize=25)

par(mar=c(5,5,4,8))

n <- 5
plot(nes[,'LNCaP--10_DAYSvLNCaP_RESIDUAL--27_DAYS'], nes[,'LNCaP--10_DAYSvLNCaP_R-CLONES--27_DAYS'],
	 col=squash::cmap(nes[,'LNCaP--10_DAYSvLNCaP_RESIDUAL--10_DAYS'], nes.cmap), pch=19,
	 xlim = c(-n,n), ylim=c(-n,n),
	 xlab='Untreated vs Residual 27 days', ylab='Untreated vs Resistant 27 days',
	 main='msVIPER Analysis of PCA Dataset', cex.lab=1.5, cex.main=1.7)

p <- 0.05
abline(h=citrus::p2z(p), lty=2)
abline(h=-citrus::p2z(p), lty=2)
abline(v=citrus::p2z(p), lty=2)
abline(v=-citrus::p2z(p), lty=2)
abline(a=0, b=1)

squash::vkey(x=n+2.2, y=-1.5, nes.cmap, 'Untreated vs.\nResidual 10 days', skip=10, stretch=0.5)

text(nes[genes,'LNCaP--10_DAYSvLNCaP_RESIDUAL--27_DAYS']+0.5, nes[genes,'LNCaP--10_DAYSvLNCaP_R-CLONES--27_DAYS'],
	 labels=genes, pos=3)

dev.off()

#############################################
########## 2.3 Transition plot
#############################################

### Input and output
rm(list=ls())
infiles <- c('f2-msviper.dir/prem_prostate-msviper_table.rda', 'f2-msviper.dir/prem_prostate-msviper_significant_gene_enrichment.rda')
outfile <- 'plots/f2-msviper/prem_prostate-msviper_transitions.png'

# Function for nice names
nicenames <- function(x) {x <- gsub('--', ' ', x); x <- gsub('_', ' ', x); x <- gsub('v', ' vs ', x); return(x)}

### Code
# Load libraries
library(citrus)

# Load data
load(infiles[1])
load(infiles[2])
nes <- -nes[complete.cases(nes),]
colnames(zscore_table) <- gsub('HALLMARK_', '', colnames(zscore_table))
zscore_table <- zscore_table[c(4, 3, 12, 11, 10, 9), c('MYC_TARGETS_V1','E2F_TARGETS','G2M_CHECKPOINT','DNA_REPAIR','APICAL_JUNCTION')]
colnames(zscore_table) <- gsub('_', ' ', colnames(zscore_table))
zscore_table <- as.matrix(t(zscore_table))

# Labels
groups <- c('LNCaP--10_DAYSvLNCaP_RESIDUAL--10_DAYS', 'LNCaP_RESIDUAL--10_DAYSvLNCaP_RESIDUAL--27_DAYS', 'LNCaP_RESIDUAL--10_DAYSvLNCaP_R-CLONES--27_DAYS')
labels <- setNames(groups, nicenames(groups))

# Count significant genes
counts <- apply(nes, 2, function(x) c(sum(x > p2z(0.05)), sum(x < -p2z(0.05))))

# Make plot
png(outfile, height=500, width=1000)

par(mar=c(4,0,4,2))

plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(-0.1,1.2), ylim=c(-0.5, 1.3))

##### 1. Black markers #####
text(x=0, y=0.5, labels='LNCaP\nUntreated', adj=0.5, cex=1.3)
text(x=0.375, y=0.5, labels='LNCaP Residual (10 days)', adj=0.3, cex=1.3)
text(x=1.02, y=0.95, labels='LNCaP Residual (27 days)', adj=0.3, cex=1.3)
text(x=1.02, y=0.05, labels='LNCaP Resistant (27 days)', adj=0.3, cex=1.3)

# Black arrows
shape::Arrows(0.1, 0.5, 0.25, 0.5, col='black', arr.type='triangle', lwd=3)
segments(0.59, 0.5, 0.63, 0.5, lwd=3)
segments(0.63, 0.5, 0.65, 0.95, lwd=3)
segments(0.63, 0.5, 0.65, 0.05, lwd=3)
shape::Arrows(0.65, 0.95, 0.9, 0.95, col='black', arr.type='triangle', lwd=3)
shape::Arrows(0.65, 0.05, 0.9, 0.05, col='black', arr.type='triangle', lwd=3)

##### 2. Colored arrows #####
# Settings
pos <- c(0.12, 0.63, 0.12, 0.75)
s <- 0.1
arrow_opts <- list('1'=c(0, 0), '2'=c(0.6, 0.45), '3'=c(0.6, -0.45))
names(arrow_opts) <- groups

for (l in names(arrow_opts)
	){
	# Get option
	opt <- arrow_opts[[l]]

	# Arrows
	shape::Arrows(pos[1]+opt[1], pos[2]+opt[2], pos[3]+opt[1], pos[4]+opt[2], col='red3', arr.type='triangle', lwd=2)
	shape::Arrows(pos[1]+opt[1]+s, pos[4]+opt[2]+0.02, pos[3]+opt[1]+s, pos[2]+opt[2]+0.02, col='navy', arr.type='triangle', lwd=2)

	# Text
	text(x=pos[1]+opt[1], y=pos[2]+opt[2]-0.07, labels=paste(counts[1,l], 'genes'))
	text(x=pos[1]+opt[1]+s, y=pos[2]+opt[2]-0.07, labels=paste(counts[2,l], 'genes'))

}

##### 3. Color heatmaps #####
# Make color map
enrich.cmap <- squash::makecmap(zscore_table, n=30, colFn=colorRampPalette(c('white','pink','red','red3')))
zscore_colors <- squash::cmap(zscore_table, enrich.cmap)

for (i in 1:length(arrow_opts))
{
	# Get option
	opt <- arrow_opts[[i]]

	# Plot colors
	n <- 2*i-1
	points(x=rep(pos[1]+opt[1], nrow(zscore_colors)), y=seq(pos[1]+opt[2]-0.1, pos[1]+opt[2]+0.23, length.out=nrow(zscore_colors)), pch=22, bg=zscore_colors[,n], cex=3.5)
	points(x=rep(pos[1]+opt[1]+s, nrow(zscore_colors)), y=seq(pos[1]+opt[2]-0.1, pos[1]+opt[2]+0.23, length.out=nrow(zscore_colors)), pch=22, bg=zscore_colors[,n+1], cex=3.5)

	# Plot text
	text(x=rep(pos[1]+opt[1], nrow(zscore_colors))+0.13, y=seq(pos[1]+opt[2]-0.1, pos[1]+opt[2]+0.23, length.out=nrow(zscore_colors)), labels=rownames(zscore_colors), pos=4)

}

# Key
squash::hkey(enrich.cmap, 'Geneset Enrichment Z-score', x=0, y=-0.45, skip=10)

# Add title
mtext('Androgen Deprivation of LNCaP - msVIPER Analysis', cex=2)

dev.off()


#############################################
########## 2.4 Zscore heatmap
#############################################

### Input and output
infile <- 'f2-msviper.dir/prem_prostate-msviper_significant_gene_enrichment.rda'
outfile <- 'plots/f2-msviper/prem_prostate-msviper_transitions.png'

### Plot
load(infile)


#######################################################
#######################################################
########## S1. Androgen Independence Analysis
#######################################################
#######################################################

#############################################
########## 1.1 AR and PSA Differential Expr.
#############################################

# Make directory
dir.create('plots/f1-androgen_independence_analysis')

### Input and Output
infiles <- c('f1-data.dir/design_table.txt', 'f1-data.dir/prem_prostate-vst.rda')
outfile <- 'plots/f1-androgen_independence_analysis/f1-ar_psa_de.png'

### Plot
# Load infiles
design_df <- tread(infiles[1])
load(infiles[2])

# Filter design dataframe by removing DMSO samples
design_df <- design_df[design_df$drug != 'DMSO',]
design_df <- design_df[design_df$group %in% c('LNCaP--10_DAYS', 'LNCaP_RESIDUAL--10_DAYS', 'LNCaP_RESIDUAL--27_DAYS', 'LNCaP_R-CLONES--27_DAYS'),]

# Get subset of KLK3 and AR expression
vsd_subset_melt <- melt(vsd_df[c('KLK3','AR','TP53'),])
colnames(vsd_subset_melt) <- c('gene_symbol', 'well', 'vst_expression')

# Replace KLK3 with PSA
vsd_subset_melt$gene_symbol <- gsub('KLK3', 'PSA', vsd_subset_melt$gene_symbol)

# Merge with design
merged_expression_df <- merge(vsd_subset_melt, design_df, by='well')

# Reorder levels
merged_expression_df$group <- factor(merged_expression_df$group, c('LNCaP--10_DAYS', 'LNCaP_RESIDUAL--10_DAYS', 'LNCaP_RESIDUAL--27_DAYS', 'LNCaP_R-CLONES--27_DAYS'))

# Plot
png(outfile, height=1500, width=600, pointsize=19)

par(mfrow=c(3,1))

for (gene in c('AR','PSA','TP53'))
{
	# Get subset for plotting
	plot_df <- merged_expression_df[merged_expression_df$gene_symbol == gene,]

	# Plot
	boxplot(vst_expression ~ group, data=plot_df, main=gene, xaxt='n', ylab='Expression Level (VST)')

	# Fix axis
	axis(1, 1:4, c('Untreated', 'Residual\n(10 days)', 'Residual\n(27 days)', 'Resistant\n(27 days)'), padj=1)
}

dev.off()

#############################################
########## 1.2 Differential Expression NES
#############################################

### Input and output
infile <- 'f5-pathway_enrichment.dir/prem_prostate-logfc_table_enrichment.rda'
outfile <- 'plots/f1-androgen_independence_analysis/f2-de_gsea.png'

### Plot
# Load libraries
library(pheatmap)

# Load infile
load(infile)

# Remove HALLMARK from rownames of table
rownames(gsea_nes_table) <- gsub('HALLMARK_', '', rownames(gsea_nes_table))

# Select hallmarks
selected_hallmarks <- c('G2M_CHECKPOINT','OXIDATIVE_PHOSPHORYLATION','PI3K_AKT_MTOR_SIGNALING','MYC_TARGETS_V1','ANDROGEN_RESPONSE')

# Get subset
gsea_subset <- -gsea_nes_table[selected_hallmarks, c('LNCaP--10_DAYSvLNCaP_RESIDUAL--10_DAYS','LNCaP_RESIDUAL--10_DAYSvLNCaP_RESIDUAL--27_DAYS','LNCaP_RESIDUAL--10_DAYSvLNCaP_R-CLONES--27_DAYS')]

rownames(gsea_subset) <- gsub('_', ' ', rownames(gsea_subset))

# Plot
png(outfile, height=500, width=300)

pheatmap(gsea_subset, color = colorRampPalette(c("navy", "white", "firebrick3"))(50), cluster_col=FALSE, gaps_col = c(1,2,3), cluster_row=FALSE)

dev.off()

#############################################
########## 1.3 MR Plot
#############################################

### Input and output
infile <- 'f2-msviper.dir/msviper_runs/LNCaP_R-CLONES--27_DAYSvLNCaP_RESIDUAL--27_DAYS.rda'
outfile <- 'plots/f1-androgen_independence_analysis/f3-mr.png'

### Plot
# Load libraries
library(viper)

# Load infile
load(infile)

### Plot
png(outfile, height=500, width=700)
plot(msviper_results, mr=c('FOXM1','TRIP13','CHEK1','AURKA','CENPI','TOP2A','UBE2C','TYMS','CENPF'))
dev.off()


#############################################
########## 1.4 AR Network Enrichment
#############################################

### Input and output
infile <- 'f5-pathway_enrichment.dir/prem_prostate-AR_network_enrichment.rda'
outfile <- 'plots/f1-androgen_independence_analysis/f4-ar_gsea.png'

### Plot
# Load infile
load(infile)

# Remove HALLMARK from rownames of table
rownames(gsea_nes_table) <- gsub('HALLMARK_', '', rownames(gsea_nes_table))

# Select hallmarks
selected_hallmarks <- c('ANDROGEN_RESPONSE', 'G2M_CHECKPOINT', 'E2F_TARGETS','GLYCOLYSIS','P53_PATHWAY')

# Get subset
gsea_subset <- gsea_nes_table[selected_hallmarks,]# c('LNCaP--10_DAYSvLNCaP_RESIDUAL--10_DAYS','LNCaP_RESIDUAL--10_DAYSvLNCaP_RESIDUAL--27_DAYS','LNCaP_RESIDUAL--10_DAYSvLNCaP_R-CLONES--27_DAYS')]

rownames(gsea_subset) <- gsub('_', ' ', rownames(gsea_subset))

# Plot
png(outfile, height=500, width=500)

pheatmap(gsea_subset[,c(1,3,2)], color = colorRampPalette(c("navy", "white", "firebrick3"))(50), cluster_col=FALSE, gaps_col = c(1,2,3), cluster_row=FALSE)

dev.off()

#############################################
########## 1.5 MYCN GSEA
#############################################

### Input and output
infiles <- c('f3-differential_expression.dir/prem_prostate-logfc_table.rda',
		     '/ifs/data/c2b2/ac_lab/dt2539/projects/project_data/f1-genelists.dir/msigdb-genesets.rda')
outfile <- 'plots/f1-androgen_independence_analysis/f5-mycn_gsea.png'

### Plot
# Load libraries
library(citrus)

# Load infile
load(infiles[1])
load(infiles[2])

# Convert rownames
rownames(logfc) <- s2e(rownames(logfc))

# Get MYC targets
myc_targets <- c(msigdb_genesets$HALLMARK_MYC_TARGETS_V1, msigdb_genesets$HALLMARK_MYC_TARGETS_V2)

# Get GSEA
gsea_results <- gsea(as.matrix(logfc)[,4], myc_targets, method='pareto')

# Plot
png(outfile, height=600, width=1000, pointsize=19)

plot_gsea(gsea_results)

dev.off()

#############################################
########## 1.6 AR GSEA
#############################################

### Input and output
infiles <- c('f3-differential_expression.dir/prem_prostate-logfc_table.rda',
		     '/ifs/data/c2b2/ac_lab/dt2539/projects/project_data/f1-genelists.dir/msigdb-genesets.rda')
outfile <- 'plots/f1-androgen_independence_analysis/f6-ar_gsea.png'

### Plot
# Load libraries
library(citrus)

# Load infile
load(infiles[1])
load(infiles[2])

# Convert rownames
rownames(logfc) <- s2e(rownames(logfc))

# Get MYC targets
myc_targets <- c(msigdb_genesets$HALLMARK_ANDROGEN_RESPONSE)

# Get GSEA
gsea_results <- gsea(as.matrix(logfc)[,4], myc_targets, method='pareto')

# Plot
png(outfile, height=600, width=1000, pointsize=19)

plot_gsea(gsea_results)

dev.off()

#############################################
########## 1.7 FOXA1 - AR correlation
#############################################

### Input and Output
infiles <- c('f1-data.dir/design_table.txt', 'f1-data.dir/prem_prostate-vst.rda')
outfile <- 'plots/f1-androgen_independence_analysis/f6-ar_foxa1_cor.png'

### Plot
# Load infiles
design_df <- tread(infiles[1])
load(infiles[2])

# Filter design dataframe by removing DMSO samples
# design_df <- design_df[design_df$cell_line != 'DU145',]# %in% c('LNCaP--10_DAYS', 'LNCaP_RESIDUAL--10_DAYS', 'LNCaP_RESIDUAL--27_DAYS', 'LNCaP_R-CLONES--27_DAYS'),]
# design_df <- design_df[design_df$group %in% c('LNCaP--10_DAYS', 'LNCaP_RESIDUAL--10_DAYS', 'LNCaP_RESIDUAL--27_DAYS', 'LNCaP_R-CLONES--27_DAYS'),]

# Get wells
wells <- sapply(unique(design_df$cell_line), function(x) design_df[design_df$cell_line==x, 'well'])

# Get expmats
expmats <- sapply(wells, function(x) vsd_df[,x])

# Plot
png(outfile, height=300, width=1000)

par(mfrow=c(1,4))

for (cond in names(expmats))
{
	x <- expmats[[cond]]['AR',]
	y <- expmats[[cond]]['FOXA1',]

	plot(x, y, xlab='AR', ylab=ifelse(cond==names(expmats)[1],'FOXA1',''), main=cond)

	abline(lm(y~x))

	print(summary(lm(y~x)))
}


dev.off()

#############################################
########## 1.8 Resistant vs Residual NES
#############################################

### Input and Output
infile <- 'f5-pathway_enrichment.dir/prem_prostate-logfc_table_enrichment.rda'
outfile <- 'plots/f1-androgen_independence_analysis/f7-de_nes.png'

### Plot
# Load libraries
library(squash)

# Load infile
load(infile)

# Remove HALLMARK from rownames of table
rownames(gsea_nes_table) <- gsub('HALLMARK_', '', rownames(gsea_nes_table))

# Select hallmarks
selected_hallmarks <- c('MYC_TARGETS_V1','G2M_CHECKPOINT','E2F_TARGETS','P53_PATHWAY','ANDROGEN_RESPONSE')

# Get NES values
nes_values <- t(gsea_nes_table[selected_hallmarks, 'LNCaP_R-CLONES--27_DAYSvLNCaP_RESIDUAL--27_DAYS'])

# Get labels
hallmark_labels <- gsub('_', ' ', selected_hallmarks)

# Get color map
nes.cmap <- makecmap(seq(-6.5, 6.5, by=0.5), colFn=colorRampPalette(c('navy','white','firebrick3')), n=30)

# Get colors
nes_colors <- cmap(nes_values, nes.cmap)

# Plot
png(outfile, height=500, width=850)

plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0, length(nes_colors)+1), ylim=c(0, 5))
points(x=1:length(nes_colors), y=rep(1, length(nes_colors)), bg=nes_colors, pch=22, cex=23, xaxt='n', yaxt='n')
text(x=1:length(nes_colors)-0.1, y=rep(2, length(nes_colors)), labels=hallmark_labels, srt=45, adj=0, cex=1.3)

hkey(nes.cmap, wh=which(nes.cmap$breaks %in% seq(-6, 6, by=3)), x=1, stretch=0.9)

dev.off()


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
command_args <- commandArgs(trailingOnly = TRUE)
str2R(command_args[1], command_args[2], command_args[3])

