
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

