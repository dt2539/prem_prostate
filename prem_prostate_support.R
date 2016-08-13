
#################################################################
#################################################################
############### Prostate Cancer Cell Line Analysis - Prem - Support
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
source('/ifs/data/c2b2/ac_lab/dt2539/projects/support/support.R')

##### 2. Other libraries #####


#######################################################
#######################################################
########## S1. Data processing
#######################################################
#######################################################

#############################################
########## 4. Get library sizes
#############################################

get_table_xy <- function(n, ncol)
{
    # Stop if not a correct division
    if ((n/ncol)%%1 != 0)
    {
        stop('n must be a multiple of ncol!')
    }

    # Give x and y arrays
    coord_df <- data.frame(x = sapply(1:n, function(x) (x-1) %% ncol) + 1,
                           y = rev(sapply((1:n)/ncol, ceiling)))

    # Return result
    return(coord_df)
}

#######################################################
#######################################################
########## S2. Specific support
#######################################################
#######################################################

#############################################
########## 2.2 (2.4) Get up/down sets
#############################################

get_signed_genesets <- function(gene_score_table)
{
    # Get up genes
    up_genes <- sapply(colnames(gene_score_table), function(x) rownames(gene_score_table)[gene_score_table[,x] > 0])
    names(up_genes) <- paste0(names(up_genes), '__up')

    # Get down genes
    down_genes <- sapply(colnames(gene_score_table), function(x) rownames(gene_score_table)[gene_score_table[,x] < 0])
    names(down_genes) <- paste0(names(down_genes), '__down')

    # Join lists
    signed_genesets <- c(up_genes, down_genes)

    # Return result
    return(signed_genesets)
}


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

