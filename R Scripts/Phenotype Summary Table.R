##################################################
## Project: Oil MP Manuscript
## Script purpose: Make a phenotype summary table to provide 
## basic statistics by location and population.
## Date: 2020-12-30
## Author: Jay Gillenwater
##################################################

# Load libraries
library(tidyverse)
library(here)
library(reshape2)
library(gt)
library(kableExtra)


# Set root directory
Rdir <- paste0(here(), "/Data/RawData")


Combined.201 <- read.csv(paste0(Rdir, "/Pop201Combined.csv"))
Combined.202 <- read.csv(paste0(Rdir, "/Pop202Combined.csv"))

Cla.201 <- read.csv(paste0(Rdir, "/Pop201Cla.csv"))
Ply.201 <- read.csv(paste0(Rdir, "/Pop201Ply.csv"))
Cas.201 <- read.csv(paste0(Rdir, "/Pop201Cas.csv"))

Cla.202 <- read.csv(paste0(Rdir, "/Pop202Cla.csv"))
Ply.202 <- read.csv(paste0(Rdir, "/Pop202Ply.csv"))
Cas.202 <- read.csv(paste0(Rdir, "/Pop202Cas.csv"))

# A function to add location and population columns to a dataframe
AddLocPop <- function(DF, Loc = "Cla", Pop = "201"){
  DF$Loc <- Loc
  DF$Pop <- Pop
  DF
}

# Add location and population columns to the phenotype data
Combined.201 <- AddLocPop(Combined.201, "Combined", "201")
Combined.202 <- AddLocPop(Combined.202, "Combined", "202")
Cla.201      <- AddLocPop(Cla.201, "Cla", "201")
Ply.201      <- AddLocPop(Ply.201, "Ply", "201")
Cas.201      <- AddLocPop(Cas.201, "Cas", "201")
Cla.202      <- AddLocPop(Cla.202, "Cla", "202")
Ply.202      <- AddLocPop(Ply.202, "Ply", "202")
Cas.202      <- AddLocPop(Cas.202, "Cas", "202")

AllData <- list(Cla.201, Ply.201, Cas.201, Cla.202, Ply.202, Cas.202, Combined.201, Combined.202)
AllData <- do.call(bind_rows, AllData) 
AllData <- AllData %>% rename(SDWT = Hundred.Seed.Weight,
                              Oil = Oil.Dry.basis,
                              Protein = Protein.Dry.basis,
                              Environment = Loc,
                              Population = Pop)

write.csv(AllData, paste0(here(), "/Data/AllPhenotypeData.csv"), row.names = FALSE)

# Split the data based on population
AllData.split <- split(AllData, AllData$Population)

# The goal of this script is to make a table to show the mean trait values for parents, and
# more detailed statistics for the phenotypic measurements across both populations
# and environments. The average values for the parents will have to be calculated first,
# then the summary information for the RILs. Following this, the two tables can then be 
# joined together

## Section: Parent data summary
##################################################

# The parents for populations 201 and 202
Parents <- c("LMN09-119", "N09-09", "LMN09-19", "N13-47")

# A function to make summary data for the parents
Parent_Summary <- function(PopulationData = AllData.split$`201`){
  PopulationData %>% 
    filter(Genotype %in% Parents) %>%
    melt(measure.vars = c("SDWT", "Oil", "Protein")) %>%
    mutate(value = round(value, 2)) %>%
    spread(Genotype, value) %>%
    rename(Trait = variable) %>%
    arrange(Population, Trait, Environment) %>%
    ungroup() %>%
    mutate_if(is.numeric, round, digits = 2)
}

# Summaries for the parents of populations 201 and 202
Parents.201.Summary <- Parent_Summary(AllData.split$`201`)
Parents.202.Summary <- Parent_Summary(AllData.split$`202`)


## Section: RIL Summary
##################################################

# A function to create simple summaries for the RILs
RIL_Summary <- function(PopulationData = AllData$`201`){
  PopulationData %>% 
    filter(!(Genotype %in% Parents)) %>%
    melt(measure.vars = c("SDWT", "Oil", "Protein")) %>%
    mutate(value = round(value, 2)) %>%
    group_by(Environment, variable) %>%
    summarise(Mean = mean(value),
              Min. = min(value),
              Max. = max(value),
              SD = sd(value)) %>%
    rename(Trait = variable) %>%
    arrange(Trait, Environment) %>%
    ungroup() %>%
    mutate_if(is.numeric, round, digits = 2)
}

RILs.201.Summary <- RIL_Summary(AllData.split$`201`)
RILs.202.Summary <- RIL_Summary(AllData.split$`202`)

## Section: Combining data
##################################################

# Join the data for each population by location and trait
SummaryData.201 <- left_join(Parents.201.Summary, RILs.201.Summary)
SummaryData.202 <- left_join(Parents.202.Summary, RILs.202.Summary)

## Section: Formatting for export
##################################################
# Remove the population column from each dataframe


SummaryGT.201 <- SummaryData.201 %>%
  select(-one_of("Population")) %>%
  gt() %>%
  tab_header(title = "Table 1: Phenotype trait values for the parental inbreds LMN09-119 and N09-09 and the derived RILs of mapping population 201 evaluated in three environments and averaged across environments. Values represent the means of parents and RILs and the minimum maximum, and standard deviation of the RILs.") %>%
  tab_spanner(label = "Parents",
              columns = vars(`LMN09-119`, `N09-09`)) %>%
  tab_spanner(label = "RILs",
              columns = vars(Mean, Min., Max., SD)) %>% 
  tab_footnote(footnote = "The environment from which the phenotype measurements come. Ply = Plymouth, Cla = Clayton. Cas = Caswell",
               locations = cells_column_labels(columns = vars(Environment))) %>% 
  tab_footnote(footnote = "The phenotype. SDWT = hundred seed weight in grams, Oil = oil content measured on a dry basis, Protein = protein content measured on a dry basis",
               locations = cells_column_labels(columns = vars(Trait))) %>%
  tab_footnote(footnote = "The parents of a population. Values are the average for each trait within each environment.",
               locations = cells_column_spanners("Parents")) %>%
  tab_footnote(footnote = "Summary statistics for the population RILs.",
               locations = cells_column_spanners("RILs"))


SummaryGT.202 <- SummaryData.202  %>%
  select(-one_of("Population")) %>% 
  gt() %>%
  tab_header(title = "Table 2: Phenotype trait values for the parental inbreds LMN09-19 and N13-47 and the derived RILs of mapping population 202 evaluated in three environments and averaged across environments. Values represent the means of parents and RILs and the minimum maximum, and standard deviation of the RILs.") %>%
  tab_spanner(label = "Parents",
              columns = vars(`LMN09-19`, `N13-47`)) %>%
  tab_spanner(label = "RILs",
              columns = vars(Mean, Min., Max., SD)) %>% 
  tab_footnote(footnote = "The environment from which the phenotype measurements come. Ply = Plymouth, Cla = Clayton. Cas = Caswell",
               locations = cells_column_labels(columns = vars(Environment))) %>% 
  tab_footnote(footnote = "The phenotype. SDWT = hundred seed weight in grams, Oil = oil content measured on a dry basis, Protein = protein content measured on a dry basis",
               locations = cells_column_labels(columns = vars(Trait))) %>%
  tab_footnote(footnote = "The parents of a population. Values are the average for each trait within each environment.",
               locations = cells_column_spanners("Parents")) %>%
  tab_footnote(footnote = "Summary statistics for the population RILs.",
               locations = cells_column_spanners("RILs"))

# Save each table as a .rtf file so that it can be added to a word document
gtsave(SummaryGT.201, paste0(here(), "\\Tables\\Pop201PhenoSummary.rtf"))
gtsave(SummaryGT.202, paste0(here(), "\\Tables\\Pop202PhenoSummary.rtf"))



  

