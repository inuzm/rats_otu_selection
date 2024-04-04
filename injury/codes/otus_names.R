#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-I", "--input_dir"), type = "character", default = "../results/", 
                help = "input directory [default= %default]", metavar = "character"),
    make_option(c("-d", "--tax_data"), type = "character", default = "https://github.com/Jonnyzamudio7/TBI/raw/main/Data/", 
                help = "path where the taxonomic dictionaries are stored [default= %default]", metavar = "character"),
    make_option(c("-O", "--out_dir"), type = "character", default = "../results/otus_taxonomic/", 
                help = "output directory [default= %default]", metavar = "character"),
    make_option(c("-m", "--model"), type = "character", default = "nb",
                help = "model to adjust: Poisson (p), Negative Binomial (nb), Zero Inflated Poisson (zip) or Zero Inflated Negative Binomial (zinb) [default= %default]",
                metavar = "integer")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# March 26th 2023
# Imanol Nuñez

#-------------------------------------------------------------------------------
# Load libraries via pacman
pacman::p_load(dplyr, tibble, tidyr)

#-------------------------------------------------------------------------------
# Prepare the data reading 
path_to_input <- opt$input_dir
path_to_output <- opt$out_dir
prefix2 <- opt$model
path_to_otus <- "significant_otus/"

# Reading the data
signif_otus_best <- read.csv(
    paste0(path_to_input, path_to_otus, prefix2, "_signif.csv")
)

# Get the different taxonomic levels present in the data
taxLevels <- unique(unlist(signif_otus_best$hlevel))

sepList <- vector("list", length = length(taxLevels))
for (i in 1:length(taxLevels)) {
    # get the taxons corresponding to a taxonomic level
    signif_temp <- signif_otus_best %>% 
        filter(hlevel == taxLevels[i])
    # Get the taxonomic names of the selected OTUs
    tempData <- read.csv(
        paste0(
            opt$tax_data, taxLevels[i], "_df_absolute.csv"
        )
    ) %>% 
        select(-c('Sample', 'Abundance')) %>% 
        filter(OTU %in% unlist(signif_temp$OTU)) %>% 
        distinct()
    # Add the names of the OTUs to the tables
    sepList[[i]] <- signif_temp %>% 
        left_join(tempData, by = "OTU")
}

# Join the data with taxonomic names
signif_otus_best_names <- sepList[[1]]
if (length(taxLevels) > 1) {
    for (i in 2:length(taxLevels)) {
        signif_otus_best_names <- signif_otus_best_names %>% 
            full_join(sepList[[i]])
    }
}

# Deleting columns that are not useful
signif_otus_best_names <- signif_otus_best_names %>% 
    select(-c(
        "pvalues", "model", "score", "adj_pvalues", "OtuClass", "g1", "g2"
    )) %>% 
    select(c('OTU', 'g1vsg2', 'hlevel', 'sign_rank', "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) %>% 
    distinct()

# Order the columns according to OTUs and the taxonomic level (useful to see if 
# nesting is present)
signif_otus_best_names %>% 
    arrange(
        OTU, 
        factor(hlevel, levels = c(
            "phylum", "order","genus"
        ))
    ) %>% 
    write.csv(
        file = paste0(
            path_to_output, "otus_", prefix2, "_tax.csv"
        ),
        row.names = FALSE
    )
# Save the produced table
