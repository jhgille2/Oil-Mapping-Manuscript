##################################################
## Project: Oil MP Manuscript
## Script purpose: Make map summary tables
## Date: 2021/01/23
## Author: Jay Gillenwater
##################################################


library(qtl)
library(tidyverse)
library(here)

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
