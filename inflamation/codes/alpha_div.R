#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-I", "--input_dir"), type = "character", default = "../results/", 
                help = "input directory [default= %default]", metavar = "character"),
    make_option(c("-d", "--tax_data"), type = "character", default = "../results/otus_taxonomic/", 
                help = "directory where the taxonomic dictionaries are stored [default= %default]", metavar = "character"),
    make_option(c("-O", "--out_dir"), type = "character", default = "../results/alpha_diversity/", 
                help = "output directory [default= %default]", metavar = "character"),
    make_option(c("-m", "--model"), type = "character", default = "nb",
                help = "model to adjust: Poisson (p), Negative Binomial (nb), Zero Inflated Poisson (zip) or Zero Inflated Negative Binomial (zinb) [default= %default]",
                metavar = "integer"),
    make_option(c("-M", "--metadata"), type = "character", default = "https://github.com/Jonnyzamudio7/TBI/raw/main/Data/metadata_sample_data.csv", 
                help = "Path to the metadata file [default = %default]", 
                metavar = "character"),
    make_option(c("-l", "--level"), type = "character", default = "order", 
                help = "hierarchical level for fill", metavar = "characer")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# March 26th 2023
# Imanol Nuñez
#-------------------------------------------------------------------------------
# Load libraries via pacman
pacman::p_load(ggplot2, ggthemes, phyloseq,                       # Plots
               plyr, reshape2, dplyr, tibble, tidyr, purrr, broom, pscl)                 # Data frame manipulation 

#-------------------------------------------------------------------------------
# Create the plot 
prefix2 <- opt$model

path_to_counts <- paste0(
    opt$input_dir, "integrated_tables/", prefix2, "_integrated.csv"
)
path_to_tax <- paste0(
    opt$input_dir, "otus_taxonomic/otus_", prefix2, "_tax.csv"
)
path_to_fig <- paste0(
    opt$out_dir, prefix2, "_relAbund.png"
)

counts <- read.csv(path_to_counts) %>% 
    select(c('ID', 'Sample', 'Abundance')) %>% 
    pivot_wider(names_from = Sample, values_from = Abundance)
taxID <- read.csv(path_to_tax)

taxData <- taxID %>% 
    mutate(ID = paste0(OTU, "_", hlevel)) %>% 
    relocate("ID") %>% 
    select(-c("OTU", "g1vsg2", "hlevel", "sign_rank")) %>% 
    group_by(ID) %>% 
    distinct() %>% 
    column_to_rownames("ID")

mdb <- read.csv(opt$metadata) %>%
    select(c('muestra', 'inflamacion')) %>%
    dplyr::rename(Sample = muestra) %>%
    dplyr::rename(Group = inflamacion) %>% 
    mutate(Group = factor(Group))

sampleData <- data.frame(
    Sample = colnames(counts)[-1]
) %>% 
    left_join(mdb, by = 'Sample') %>% 
    mutate(rowName = Sample) %>% 
    column_to_rownames('rowName')

idx <- 1:nrow(sampleData)

sampleData <- sample_data(
    sampleData[idx, ]
)

countsData <- counts[, c(1, idx + 1)] %>% 
    column_to_rownames("ID") %>% 
    otu_table(taxa_are_rows = TRUE)

phylo_sv <- phyloseq(
    countsData,
    tax_table(as.matrix(taxData)),
    sampleData
)

alphadiv <- estimate_richness(
    phylo_sv,
    measures=c("Observed", "Shannon", "Chao1")
) %>% 
    as.data.frame() %>% 
    mutate(Sample = rownames(.)) %>% 
    pivot_longer(-Sample, names_to = "index", values_to = "Value") %>% 
    mutate(Sample = gsub('[.]', '-', Sample)) %>% 
    left_join(mdb, by = 'Sample') %>% 
    filter(substr(index, 1, 2) != "se") %>% 
    ggplot(aes(x = Group, y = Value, color = Group, group = Group, fill = Group)) +
    geom_violin(alpha = 0.1) +
    geom_point(position = position_jitter(w = 0.1, h = 0)) +
    facet_wrap( index ~ ., scales = "free", ncol = 1) +
    theme_fivethirtyeight() +
    theme(axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.margin = margin()) +
    ggtitle("Alpha diversity indexes for samples with selected variables") +
    guides(colour = guide_legend(nrow = 3))

ggsave(
    filename = path_to_fig,
    plot = alphadiv,
    dpi = 180, width = 9, height = 16
)
