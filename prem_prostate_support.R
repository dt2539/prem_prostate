
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
########## 2.1 (2.4) Get enrichment
#############################################

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


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

