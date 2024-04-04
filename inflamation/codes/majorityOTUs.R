#!/usr/bin/env Rscript
#install.packages("optparse", repos = "http://cran.us.r-project.org", lib = "./libs")
#install.packages("pacman", repos = "http://cran.us.r-project.org", lib = "./libs")
#library("optparse", lib.loc = "./libs")
#library("pacman", lib.loc = "./libs")
pacman::p_load(optparse)

option_list = list(
    make_option(c("-o", "--original"), type = "character", default = TRUE, 
                help = "input directory [default= %default]", metavar = "character"),
    make_option(c("-O", "--out_dir"), type = "character", default = "../results/majorityOTUs/", 
                help = "output directory [default= %default]", metavar = "character"),
    make_option(c("-m", "--model"), type = "character", default = "nb",
                help = "adjusted model: Poisson (p), Negative Binomial (nb), Zero Inflated Poisson (zip) or Zero Inflated Negative Binomial (zinb), best [default= %default]",
                metavar = "integer"),
    make_option(c("-p", "--perc"), type = "numeric", default = 100, 
                help = "percentage of samples where OTU is located",
                metavar = "numeric")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# March 26th 2023
# Imanol Nuñez

#-------------------------------------------------------------------------------
# Load libraries via pacman
pacman::p_load(ggplot2, ggthemes,                       # Plots
               dplyr, tibble, tidyr, purrr, broom, pscl)                 # Data frame manipulation 
#install.packages("ggplot2", repos = "http://cran.us.r-project.org", lib = "./libs")
#library(ggplot2, lib.loc = "./libs")
#install.packages("ggthemes", repos = "http://cran.us.r-project.org", lib = "./libs")
#library(ggthemes, lib.loc = "./libs")
#install.packages("dplyr", repos = "http://cran.us.r-project.org", lib = "./libs")
#library(dplyr, lib.loc = "./libs")
#install.packages("tibble", repos = "http://cran.us.r-project.org", lib = "./libs")
#library(tibble, lib.loc = "./libs")
#install.packages("tidyr", repos = "http://cran.us.r-project.org", lib = "./libs")
#library(tidyr, lib.loc = "./libs")
#install.packages("purrr", repos = "http://cran.us.r-project.org", lib = "./libs")
#library(purrr, lib.loc = "./libs")
#install.packages("broom", repos = "http://cran.us.r-project.org", lib = "./libs")
#library(broom, lib.loc = "./libs")
#install.packages("pscl", repos = "http://cran.us.r-project.org", lib = "./libs")
#library(pscl, lib.loc = "./libs")

#-------------------------------------------------------------------------------
OTUs_majority_samples <- function(perc, pre_path) {
    tlevels <- c(
        "genus", "order", "phylum"
    )
    for (i in 1:length(tlevels)) {
        db1 <- read.csv(paste0(pre_path, tlevels[i], "_df_absolute.csv")) %>% 
            mutate(OTU = paste0(OTU, "_", tlevels[i])) %>%
            dplyr::rename(ID = OTU)
        if (i == 1) {
            db_all <- db1 %>% 
                filter(!is.character(ID))
        }
        majorityOTUs <- db1 %>% 
            group_by(ID) %>% 
            summarise(perc_app = mean(Abundance > 0)) %>% 
            filter(perc_app >= perc / 100)
        db1_selected <- db1 %>% 
            filter(ID %in% unlist(majorityOTUs$ID)) 
            
        db_all <- db_all %>% 
            full_join(db1_selected)
    }
    orderAux <- db_all %>% 
        group_by(ID) %>% 
        summarise(totAb = sum(Abundance)) %>% 
        arrange(desc(totAb))
    db_all <- db_all %>% 
        arrange(factor(ID, levels = unlist(orderAux$ID)))
    return( db_all )
}

OTUs_majority_samples_selected <- function(perc, pre_path) {
    db1 <- read.csv(paste0(pre_path, "integrated.csv"))
    majorityOTUs <- db1 %>% 
        group_by(ID) %>% 
        summarise(perc_app = mean(Abundance > 0)) %>% 
        filter(perc_app >= perc / 100)
    db1_selected <- db1 %>% 
        filter(ID %in% unlist(majorityOTUs$ID)) 
    orderAux <- db1_selected %>% 
        group_by(ID) %>% 
        summarise(totAb = sum(Abundance)) %>% 
        arrange(desc(totAb))
    db1_selected <- db1_selected %>% 
        arrange(factor(ID, levels = unlist(orderAux$ID)))
    return( db1_selected )
}

#-------------------------------------------------------------------------------
initial_path <- ifelse(
    opt$original, 
    "https://github.com/Jonnyzamudio7/TBI/raw/main/Data/",
    "../results/integrated_tables/"
)

prefix2 <- ifelse(
    opt$original,
    "",
    opt$model
)
prefix3 <- ifelse(
    opt$original,
    "",
    "integrated"
)

path_to_counts <- paste0(
    initial_path, 
    ifelse(opt$original, "", paste0(prefix2, "_"))
)

if (opt$original) {
    majorOTUs <- OTUs_majority_samples(
        opt$perc, path_to_counts
    )
} else {
    majorOTUs <- OTUs_majority_samples_selected(
        opt$perc, path_to_counts
    )
}

path_to_out_counts <- paste0(
    opt$out_dir, 
    ifelse(
        opt$original, 
        "original", 
        paste0("selected_", opt$model)
    ), 
    "_", opt$perc, ".csv"
)

write.csv(x = majorOTUs, file = path_to_out_counts, row.names = FALSE)