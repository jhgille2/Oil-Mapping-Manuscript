##################################################
## Project: Oil MP Manuscript
## Script purpose: Create an ANOVA table for all oil MP data that can be used in the manuscript
## Date: 2021-01-17
## Author: Jay Gillenwater
##################################################

# Load libraries
library(tidyverse)
library(here)
library(broom)
library(gt)

# Read in the phenotype data
AllPheno <- read_csv(paste0(here(), "/Data/RawData/FullData.csv")) %>%
  group_by(Population) %>%
  dplyr::filter(Rep != 3, Genotype != "N18-1846") %>%
  mutate(Genotype = factor(Genotype),
         Loc = factor(Loc),
         Rep = factor(Rep)) %>%
  ungroup()

# A general linear model to map to each data set
PhenoFn <- function(df){
  lm(value ~ Genotype*Loc + Rep:Loc, data = df)
}

# A function to produce a tidy analysis of variance
Tidyaov <- function(Model){
  tidy(aov(Model))
}

# Fit models, get analysis of variance, and format to a presentation-ready table
AllPheno_models <- AllPheno %>% 
  pivot_longer(cols = c(8:10)) %>%        # Pivot the data from the phenotypes
  group_by(Population, name) %>%          # Group the phenotype data by population and trait. Models will be fit within each of these groups
  nest() %>%                              # Nest the data using these groups
  mutate(model = map(data, PhenoFn),      # Map the modeling function, and the anova tidying function
         AOV   = map(model, Tidyaov),     
         name  = recode(name,             # Recode the traits            
                        Oil.Dry.basis       = "Oil", 
                        Hundred.Seed.Weight = "SDWT", 
                        Protein.Dry.basis   = "Protein")) %>%
  dplyr::select(-one_of(c("data", "model"))) %>%
  rename(Trait = name) %>%                # Rename the "name" column
  unnest(AOV) %>%                         # Unnest the analysis of variance data
  ungroup()


AllPheno %>% 
  pivot_longer(cols = c(8:10)) %>%        # Pivot the data from the phenotypes
  group_by(Population, name) %>%
  count(name)
