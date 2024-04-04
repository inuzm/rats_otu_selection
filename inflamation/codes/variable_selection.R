#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-I", "--input_dir"), type = "character", default = "https://github.com/Jonnyzamudio7/TBI/raw/main/Data/", 
                help = "input directory [default= %default]", metavar = "character"),
    make_option(c("-O", "--out_dir"), type = "character", default = "../results/", 
                help = "output directory [default= %default]", metavar = "character"),
    make_option(c("-l", "--level"), type = "character", default = "genus", 
                help = "On which hierarchical level will the selection run (genus, order, phylum) [default= %default]",
                metavar = "character"),
    make_option(c("-k", "--best_k"), type = "integer", default = 5, 
                help = "number of significant OTUs [default= %default]", metavar = "integer"),
    make_option(c("-m", "--model"), type = "character", default = "nb",
                help = "model to adjust: Poisson (p), Negative Binomial (nb), Zero Inflated Poisson (zip) or Zero Inflated Negative Binomial (zinb) [default= %default]",
                metavar = "integer")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# March 25th, 2024
# Imanol Nuñez

#-------------------------------------------------------------------------------
# Load libraries via pacman
pacman::p_load(ggplot2, ggthemes,                       # Plots
               dplyr, tibble, tidyr, purrr, broom, pscl)                 # Data frame manipulation 
# library(ggplot2)
# library(ggthemes)
# library(dplyr)
# library(tibble)
# library(tidyr)
# library(purrr)
# library(broom)
# library(pscl)

#-------------------------------------------------------------------------------
# Safety functions to keep differentialOtusPvalues running even if certain 
# fits cannot be achieved by numerical errors
poss_glm <- possibly(.f = glm, otherwise = NULL)
poss_glm.nb <- possibly(.f = MASS::glm.nb, otherwise = NULL)
poss_zinfl <- possibly(.f = pscl::zeroinfl, otherwise = NULL)
poss_AIC <- possibly(.f = AIC, otherwise = Inf)
poss_BIC <- possibly(.f = BIC, otherwise = Inf)

#-------------------------------------------------------------------------------
# Compute the AIC score for a list of models
modelCountsScore <- function(models, method = "AIC") {
    if (method == "AIC") {
        scores <- poss_AIC(models)
    } else {
        scores <- poss_BIC(models)
    } 
    return(scores)
}

#-------------------------------------------------------------------------------
# Extract the p-value for a non-constant term in a regression model of the 
# form y ~ 1 + x if it is not zero inflated, and a regression of the form 
# y ~ (1 + x | x) if a model is zero inflated
pValueFromSummary <- function(fit, zi = FALSE) {
    if (is.null(fit)) {
        pValue <- NA
    } else {
        if (zi) {
            pValue <- summary(fit)$coefficients$count[,4][2]
        } else {
            pValue <- summary(fit)$coefficients[,4][2]
        }
    }
    return(pValue)
}

#-------------------------------------------------------------------------------
# Select the best model for the counts of an OTU
# and return the associated p-value and name of the selected model
modelFitting <- function(db, formula, model, method = "AIC") {
    # Fit Poisson, Negative Binomial, Zero Inflated Poisson and 
    # Zero Inflated Negative Binomial models
    if (model == "p") {
        modelFit <- poss_glm(
            formula = formula, 
            family = "poisson", 
            data = db
        )
        zi_mod <- FALSE
    } else if (model == "nb") {
        modelFit <- poss_glm.nb(
            formula = formula, 
            data = db
        )
        zi_mod <- FALSE
    } else if (model == "zip") {
        modelFit <- poss_zinfl(
            formula = formula, 
            dist = "poisson", 
            data = db
        )
        zi_mod <- TRUE
    } else if (model == "zinb") {
        modelFit <- poss_zinfl(
            formula = formula, 
            dist = "negbin", 
            data = db
        )
        zi_mod <- TRUE
    }
    # Compute the scores according to a method (AIC by default)
    scores <- modelCountsScore(models = modelFit, method = method)
    # Compute the p-value
    pValueModel <- pValueFromSummary(
        fit = modelFit, 
        zi = zi_mod
    )
    return(
        list(pvalues = pValueModel, model = model, score = scores)
    )
}

#-------------------------------------------------------------------------------
# differentialOtusPvalues calculates the p-values according to log2-fold change 
# for every pair of rat groups
# This allows us to decide if an OTU is differentially abundant
# This only works correctly for absolute frequence 

differentialOtusPvalues <- function(db, model) {
    #####
    # Get pvalues fitting a count model to db
    
    ### First approximation
    # We assume that the data is a table with three columns: the OTUS, the 
    # sample name and the number of reads of an OTU in an sample
    # Then we add a column with the number of reads for each sample and a 
    # logarithmic transformation
    # We finally add the city and year from where the sample is from
    db <- db %>% 
        as_tibble() %>%  
        group_by(Sample) %>% 
        mutate(nReads = sum(Abundance)) %>% 
        mutate(
            logNreads = ifelse(
                nReads == 0, 1e-30, log(nReads)
            ))
    
    # Get all the groups in the sample
    groups_rats <- unique(unlist(db$Group))
    ngroups <- length(groups_rats)
    
    # Initialize the data frame where p-values, and adjusted p-values, 
    # will be stored
    # As we will be adjusting a model for each OTU, and each pair of 
    # city-year, this IDs will also be stored
    pValues <- data.frame(
        OTU = character(), pvalues = double(), model = character(), 
        score = double(), g1 = integer(), 
        g2 = integer(), g1vsg2 = character(), adj_pvalues = double()
    )
    pb <- txtProgressBar(min = 0, max = ngroups * (ngroups - 1) / 2, style = 3)
    k <- 1
    # Construct the formula for the models
    if (model %in% c("zip", "zinb")) {
        formulaModels <- formula(Abundance ~ offset(logNreads) + Group | Group)
    } else {
        formulaModels <- formula(Abundance ~ offset(logNreads) + Group)
    }
    for (i in 1:(ngroups - 1)) {
        for (j in (i+1):ngroups) {
            # Given two rats groups, we first filter the data that 
            # corresponds to said samples, which we convert to a factor with 
            # levels 0 and 1
            # We then eliminate any OTU that had 0 reads among all the sample 
            # The nesting is done at the OTU level, to adjust a negative 
            # binomial model for each OTU
            # After the models are adjusted, we extract the p-value of the 
            # coefficient corresponding to the dummy variable given by the 
            # rat group
            # After unnesting the p-values, we only keep the OTUs ID and their 
            # p-values, adding for each group, a column with the ID of 
            # group being compared
            # We then filter any OTUs such that the computed p-value is NA, 
            # to then compute the adjusted p-values via the Benjamini-Hochberg
            # correction
            tempPvalues <- db %>% 
                filter(Group %in% groups_rats[c(i, j)]) %>% 
                mutate(Group = factor(Group)) %>% 
                ungroup() %>% group_by(OTU) %>% 
                mutate(otuReads = sum(Abundance)) %>% 
                ungroup() %>% 
                filter(otuReads > 0) %>% 
                select(-otuReads) %>% 
                nest(-OTU) %>%  
                mutate(
                    modelFit = map(data, ~modelFitting(
                        formula = formulaModels, 
                        model = model, 
                        db = .
                    ))
                ) %>% 
                unnest_wider(modelFit) %>% 
                select(c("OTU", "pvalues", "model", "score")) %>% 
                mutate(
                    g1 = groups_rats[i], 
                    g2 = groups_rats[j], 
                    g1vsg2 = paste0(g1, "_vs_", g2)
                ) %>% 
                filter(!is.na(pvalues)) %>% 
                filter(is.finite(score)) %>% 
                mutate(
                    adj_pvalues = p.adjust(pvalues, method = "fdr")
                ) 
            # We join the computed p-values to the data.frame that contains 
            # the computed p-values in previous iterations
            pValues <- pValues %>% 
                full_join(tempPvalues, by = c(
                    "OTU", "pvalues", "model", "score", "g1", "g2", "g1vsg2", "adj_pvalues"
                ))
            setTxtProgressBar(pb, k)
            k <- k + 1
        }
    }
    
    return(pValues)
}

#-------------------------------------------------------------------------------
# computePvaluesLevel computes a matrix of p-values for the negative binomial 
# model adjusted by differentialOtusPvalues for the count data of 
# reads for the hierarchical level hLevel
# This function may take into account a set of indices for training 
computePvaluesLevel <- function(hLevel, path_to_counts, 
                                path_to_pvalues = NULL, model = "nb") {
    # Test if the count data is for reads or for assembly data
    db <- read.csv(paste0(path_to_counts, hLevel, "_df_absolute.csv")) %>% 
        select(c('OTU', 'Sample', 'Abundance', 'inflamacion')) %>% 
        rename(Group = inflamacion)
    ### It can also be select(-c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'))
    ### but as a first approximation we will consider a 'simpler' model
    # compute p-values
    tempPvalues <- differentialOtusPvalues(db, model)
    # Save p-values for further exploratory analysis
    #write.csv(
    #    tempPvalues, 
    #    file = paste0(path_to_pvalues, "train_pvalues_", hLevel, ".csv"), 
    #    row.names = FALSE
    #)
    # Add the level of the p_value
    tempPvalues <- tempPvalues %>% mutate(hlevel = hLevel)
    return(tempPvalues)
}

#-------------------------------------------------------------------------------
# Given a table of p-values, getKOtus gets the k most significant OTUs 
# for each pair of city/year, specifying the hierarchical level of said OTU
getKOtus <- function(db, k) {
    # Put all p-values into a single column, identifying them by the cities 
    # being contrasted
    #    pValueslong <- db %>% 
    #    pivot_longer(cols = -c("OTU", "hlevel"), 
    #                 names_to = "cities", 
    #                 values_to = "p-value") %>% 
    #    mutate(city1 = substr(cities, 1, 5), city2 = substr(cities, 10, 15))
    
    # Which comparisons were made
    comparisons <- unique(unlist(db$g1vsg2))
    # Initialize the data frame for the k most significant OTUs per 
    # city vs city contrast
    reducedK <- data.frame(
        OTU = character(), pvalues = double(), model = character(), 
        score = double(), g1 = integer(), 
        g2 = integer(), g1vsg2 = character(), adj_pvalues = double(),  
        hlevel = character(), sign_rank = integer()
    )
    # Get the significant OTus
    for(i in 1:length(comparisons)) {
        tempData <- db %>% 
            filter(g1vsg2 == comparisons[i]) %>% 
            arrange(adj_pvalues) %>% 
            head(n = k) %>% 
            mutate(sign_rank = 1:k)
        reducedK <- reducedK %>% 
            full_join(tempData, by = c(
                "OTU", "pvalues", "model", "score", "g1", "g2", "g1vsg2", "adj_pvalues",
                "hlevel", "sign_rank"
            ))
    }
    
    # Add an ID for OTU - class
    reducedK <- reducedK %>% 
        mutate(OtuClass = paste(OTU, hlevel, sep = ""))
    
    return(reducedK)
}

#-------------------------------------------------------------------------------
# Given a list of significant OTUs, construct the integrated sample with only 
# these OTUs. 
constructIntegratedData <- function(sign_otus, path_to_counts) {
    hLevels <- unique(sign_otus$hlevel)
    # Initialize a list where the subsetting will be done by taxonomical levels
    dataSubsets <- vector("list", length(hLevels))
    for (i in 1:length(hLevels)) {
        tempData <- read.csv(
            paste0(path_to_counts, hLevels[i], "_df_absolute.csv")
        )
        # Subset which significant OTUs correspond to a given taxonomical 
        # level
        tempOTUs <- sign_otus %>%  
            filter(hlevel == hLevels[i])
        # Subset those significant OTUs from the count data 
        dataSubsets[[i]] <- tempData[tempData[, 1] %in% tempOTUs$OTU, ] %>% 
            mutate(hlevel = hLevels[i])
    }
    # Merge all of the subsetted data
    retDF <- dataSubsets[[1]]
    if (length(hLevels) > 1) {
        for (i in 2:length(hLevels)) {
            retDF <- retDF %>% 
                full_join(dataSubsets[[i]])
        }
    }
    # 
    retDF <- retDF %>% 
        mutate(ID = paste0(OTU,"_", hlevel)) %>% 
        dplyr::select(-c("OTU", "hlevel")) %>% 
        relocate(ID)
    return(retDF)
}

#-------------------------------------------------------------------------------
# The function variableSelection implements all of the steps to select a subset 
# of OTUs that allow us to differentiate between cities
variableSelection <- function(hlevels = c("phylum", "order", "genus"), 
                              path_to_counts, kpvalues = 5, 
                              model = "nb") {
    # Initialize a list to save the matrices of p-values for every run
    pValuesList <- vector("list", length = length(hlevels))
    for (i in 1:length(hlevels)) {
        cat(sprintf("Starting with %s\n", hlevels[i]))
        pValuesList[[i]] <- computePvaluesLevel(
            hLevel = hlevels[i], 
            path_to_counts = path_to_counts, 
            model = model
        )
        cat(sprintf("\n%d of %d done\n", i, length(hlevels)))
    }
    # Merge all of the p-values matrices
    pValues <- pValuesList[[1]]
    if (length(hlevels) > 1) {
        for (i in 2:length(hlevels)) {
            pValues <- pValues %>% 
                full_join(pValuesList[[i]])
        }
    }
    # Given the complete list of p-values, get the k most significant for every 
    # pair of cities
    significantOtus <- getKOtus(pValues, kpvalues)
    # Construct the integrated data with the significant OTUs 
    integratedTable <- constructIntegratedData(significantOtus,
                                               path_to_counts)
    return(list(pValues, significantOtus, integratedTable))
}

#-------------------------------------------------------------------------------
# Test
path_to_counts <- opt$input_dir
kpvalues <- opt$best_k
hlevels <- c("phylum", "order", "genus")


smth <- variableSelection(path_to_counts = path_to_counts,
                          kpvalues = kpvalues, 
                          hlevels = hlevels, model = opt$model)

write.csv(
    smth[[1]],  
    file = paste0(
        opt$out_dir, "pValues/", opt$model, "_", "pvalues", ".csv"
    ),
    row.names = FALSE
)
write.csv(
    smth[[2]],  
    file = paste0(
        opt$out_dir, "significant_otus/", opt$model, "_", "signif", ".csv"
    ),
    row.names = FALSE
)
write.csv(
    smth[[3]],  
    file = paste0(
        opt$out_dir, "integrated_tables/", opt$model, "_", "integrated", ".csv"
    ),
    row.names = FALSE
)

pplot1 <- smth[[1]] %>% 
    ggplot(aes(x = g1vsg2, y = -log(adj_pvalues), colour = hlevel)) + 
    geom_hline(yintercept = -log(1e-3), colour = "hotpink") + 
    geom_point(alpha = 0.5, size = 1) + 
    theme_few() + 
    ylab("-log(p-value)") + 
    xlab("") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(0.5)),
          legend.position = "top") + 
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4)))

ggsave(
    plot = pplot1, 
    filename = paste0(
        opt$out_dir, "pValues/", opt$model, "_", "log_pvalues", ".png"
    ),
    dpi = 180, width = 12, height = 6.75
)
