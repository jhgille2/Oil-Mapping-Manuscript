##################################################
## Project: Oil MP Manuscript
## Script purpose: Make map summary tables
## Date: 2021/01/23
## Author: Jay Gillenwater
##################################################


library(qtl)
library(tidyverse)
library(here)
library(gt)

# Chromosome to LG conversion
ChrtoLG <- read.csv("https://raw.githubusercontent.com/jhgille2/SoybaseData/master/SoybaseLGAssignments.csv", stringsAsFactors = FALSE)

# A function to make a summary table of a linkage map by chromosome
MapSummary <- function(Cross){
  
  # Get the map as a table
  MapTable <- pull.map(Cross, as.table = TRUE)
  MapTable$Marker <- rownames(MapTable) # Make a column for the marker names
  
  
  # A function to format LG names, converts chromosome numbers to a 
  # string with the chromosome number and LG name
  ChrNameFormat <- function(ChrNum){
    LG <- ChrtoLG$LG[which(ChrtoLG$Chr == ChrNum)]
    
    paste(ChrNum, " (", LG, ")", sep = "")
  }
  
  MapSummaryTable <- MapTable %>% 
    group_by(chr) %>%
    summarise(nMarks = n(),
              MapLen = max(pos),
              AverageSpacing = mean(diff(pos)))
  
  MapSummaryTable$Chromosome <- as.character(ChrNameFormat(as.numeric(as.character(MapSummaryTable$chr))))
  
  # Calculate average and total for marker numbers and lg length
  SummaryFunc_avg <- function(x) if(is.numeric(x)) mean(x) else ''
  SummaryFunc_sum <- function(x) if(is.numeric(x)) sum(x) else ''
  
  MapSummaryTable$AverageSpacing <- round(MapSummaryTable$AverageSpacing, 2)
  
  
  MapSummaryTable <- MapSummaryTable %>% dplyr::select(Chromosome, nMarks, MapLen, AverageSpacing)
  MapSummaryTable <- as.data.frame(MapSummaryTable)
  
  AvgRow                         <- as.data.frame(lapply(MapSummaryTable, SummaryFunc_avg))
  AvgRow$AverageSpacing          <- round(AvgRow$AverageSpacing, 2)
  MapSummaryTable$AverageSpacing <- as.character(MapSummaryTable$AverageSpacing)
  AvgRow$AverageSpacing          <- as.character(AvgRow$AverageSpacing)
  
  SumRow   <- as.data.frame(lapply(MapSummaryTable, SummaryFunc_sum))
  
  AvgRow[, 1] <- as.character(AvgRow[, 1])
  SumRow[, 1] <- as.character(SumRow[, 1])
  
  AvgRow[1, 1] <- 'Mean'
  SumRow[1, 1] <- 'Total'
  
  MapSummaryTable[nrow(MapSummaryTable) + 1,] <- NA
  MapSummaryTable <- bind_rows(MapSummaryTable, bind_rows(AvgRow, SumRow))
  
  
  MapSummaryTable$MapLen         <- round(MapSummaryTable$MapLen, 2)
  MapSummaryTable$nMarks         <- round(MapSummaryTable$nMarks, 1)
  
  
  colnames(MapSummaryTable) <- c("Chromosome",
                                 "Number of Markers",
                                 "Chromosome size (cM)",
                                 "Average marker spacing")
  
  MapSummaryTable[is.na(MapSummaryTable)] <- ''
  MapSummaryTable
}

# Read in crosses
load(paste0(here(), "/Data/AllCrosses.RData"))

# Summary tables for each map
Pop201MapSummary <- MapSummary(AllCrosses[[1]])
Pop202MapSummary <- MapSummary(AllCrosses[[8]])


# Formatting and export
SummaryMapGT.201 <- Pop201MapSummary %>%
  gt() %>%
  tab_header(title = "Table S1: Summary data for the linkage map created for soybean oil mapping population 201 using 421 SNPs genotyped using the Illumina Infinium SoySNP6K BeadChip.") %>%
  tab_footnote(footnote = "Chromosme number and the corresponding linkage group name for each linkage group in the linkage map.",
               locations = cells_column_labels(columns = vars(Chromosome))) %>% 
  tab_footnote(footnote = "The number of SNP markers which were positioned on each linkage group.",
               locations = cells_column_labels(columns = vars(`Number of Markers`))) %>%
  tab_footnote(footnote = "The size of each linkage group in centimorgans. Genetic distances were calculated using Kosambi's function.",
               locations = cells_column_labels(columns = vars(`Chromosome size (cM)`))) %>%
  tab_footnote(footnote = "The average spacing betweebneach marker of a linkage group in centimorgans.",
               locations = cells_column_labels(columns = vars(`Average marker spacing`)))

SummaryMapGT.202 <- Pop202MapSummary %>%
  gt() %>%
  tab_header(title = "Table S2: Summary data for the linkage map created for soybean oil mapping population 202 using 416 SNPs genotyped using the Illumina Infinium SoySNP6K BeadChip.") %>%
  tab_footnote(footnote = "Chromosme number and the corresponding linkage group name for each linkage group in the linkage map.",
               locations = cells_column_labels(columns = vars(Chromosome))) %>% 
  tab_footnote(footnote = "The number of SNP markers which were positioned on each linkage group.",
               locations = cells_column_labels(columns = vars(`Number of Markers`))) %>%
  tab_footnote(footnote = "The size of each linkage group in centimorgans. Genetic distances were calculated using Kosambi's function.",
               locations = cells_column_labels(columns = vars(`Chromosome size (cM)`))) %>%
  tab_footnote(footnote = "The average spacing between each marker of a linkage group in centimorgans.",
               locations = cells_column_labels(columns = vars(`Average marker spacing`)))


# Save each table as a .rtf file so that it can be added to a word document
gtsave(SummaryMapGT.201, paste0(here(), "\\Tables\\Pop201LinkageMapSummary.rtf"))
gtsave(SummaryMapGT.202, paste0(here(), "\\Tables\\Pop202LinkageMapSummary.rtf"))
