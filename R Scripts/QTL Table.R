##################################################
## Project: Oil MP Manuscript
## Script purpose: Reformat the QTL table for more concise presentation
## Date: 2020-12-31
## Author: Jay Gillenwater
##################################################

# Load libraries
library(tidyverse)
library(reshape2)
library(here)
library(reshape2)
library(gt)
library(kableExtra)
library(readxl)
library(skimr)


# Read the full QTL Table
AllQTL <- read.csv(paste0(here(), "/Data/AllQTL.csv"))

# Fix the Consensus positions
AllQTL$ConsensusL_New <- NA
AllQTL$ConsensusR_New <- NA
AllQTL$ConsensusLPos_New <- NA
AllQTL$ConsensusRPos_New <- NA

for(i in 1:nrow(AllQTL)){
  if(AllQTL$ConsensusLPos[[i]] > AllQTL$ConsensusRPos[[i]]){
    
    AllQTL$ConsensusL_New[[i]] <- AllQTL$ConsensusR[[i]]
    AllQTL$ConsensusR_New[[i]] <- AllQTL$ConsensusL[[i]]
    
    AllQTL$ConsensusLPos_New[[i]] <- AllQTL$ConsensusRPos[[i]]
    AllQTL$ConsensusRPos_New[[i]] <- AllQTL$ConsensusLPos[[i]]
  }else{
    
    AllQTL$ConsensusL_New[[i]] <- AllQTL$ConsensusL[[i]]
    AllQTL$ConsensusR_New[[i]] <- AllQTL$ConsensusR[[i]]
    
    AllQTL$ConsensusLPos_New[[i]] <- AllQTL$ConsensusLPos[[i]]
    AllQTL$ConsensusRPos_New[[i]] <- AllQTL$ConsensusRPos[[i]]
  }
}

AllQTL$ConsensusL <- AllQTL$ConsensusL_New
AllQTL$ConsensusR <- AllQTL$ConsensusR_New

AllQTL$ConsensusLPos <- AllQTL$ConsensusLPos_New
AllQTL$ConsensusRPos <- AllQTL$ConsensusRPos_New

AllQTL_table <- AllQTL %>% 
  select(-one_of("ConsensusL_New", "ConsensusR_New", "ConsensusLPos_New", "ConsensusRPos_New")) %>%
  mutate(`Marker Interval` = paste(FlankL, FlankR, sep = "-"), `Consensus Interval` = paste(round(ConsensusLPos, 2), round(ConsensusRPos, 2), sep = "-")) %>%
  select(-one_of("FlankL", "FlankR", "ConsensusL", 'ConsensusR')) %>%
  select(Label, Loc, Pop, Trait, QTLChr, `Marker Interval`, `Consensus Interval`, LOD, var, QTLest) %>%
  rename(Environment = Loc, Population = Pop, Chromosome = QTLChr, QTL = Label, PVE = var, `Additive Effect` = QTLest)

AllQTL_table$LOD <- round(AllQTL_table$LOD, 2)

AllQTL_table_gt <- AllQTL_table  %>%
  gt() %>%
  tab_footnote(footnote = "The environment from which the phenotype measurements come. Ply = Plymouth, Cla = Clayton. Cas = Caswell.",
               locations = cells_column_labels(columns = vars(Environment))) %>% 
  tab_footnote(footnote = "The phenotype. SDWT = hundred seed weight, Oil = oil content measured on a dry basis, Protein = protein content measured on a dry basis.",
               locations = cells_column_labels(columns = vars(Trait))) %>% 
  tab_footnote(footnote = "The population which the QTL was found in.",
               locations = cells_column_labels(columns = vars(Population))) %>% 
  tab_footnote(footnote = "The chromosome which the QTL was found on.",
               locations = cells_column_labels(columns = vars(Chromosome))) %>% 
  tab_footnote(footnote = "The logarithm of odds for the QTL.",
               locations = cells_column_labels(columns = vars(LOD)))  %>% 
  tab_footnote(footnote = "Percent of the phenotypic variation explained by the QTL.",
               locations = cells_column_labels(columns = vars(PVE))) %>% 
  tab_footnote(footnote = "The additive effect of the QTL.",
               locations = cells_column_labels(columns = vars(`Additive Effect`))) %>% 
  tab_footnote(footnote = "The flanking markers of the QTL.",
               locations = cells_column_labels(columns = vars(`Marker Interval`))) %>% 
  tab_footnote(footnote = "The interval on the GmConsensus 4.0 map corresponding to the QTL flanking markers.",
               locations = cells_column_labels(columns = vars(`Consensus Interval`)))

# Save each table as a .rtf file so that it can be added to a word document
gtsave(AllQTL_table_gt, paste0(here(), "\\Tables\\QTLSummary.rtf"))


## Section: Summary information
##################################################

# The tables in this section are primarily used to write descriptions in the manuscript
# but are not included as-is

# Detailed descriptive statistics of variables
AllQTL_table %>% group_by(Trait) %>% skim()
