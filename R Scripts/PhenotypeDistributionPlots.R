##################################################
## Project: Oil MP Manuscript
## Script purpose: Functions for plotting distributions
## of phenotype(s) with the values of population parents shown
## Date: 2020-01-09
## Author: Jay Gillenwater
##################################################


## Section: Environment/data loading
##################################################

# Load packages
library(tidyverse)
library(here)
library(hrbrthemes)
library(cowplot)

# Read in the phenotype data
AllPheno <- read_csv(paste0(here(), "/Data/AllPhenotypeData.csv"))

# Filter to just the combined environment
AllPheno <- AllPheno %>%
  filter(Environment == "Combined")

# The names of the parent genotypes 
ParentGenos        <- c("LMN09-119", "N09-09", "LMN09-19", "N13-47")
names(ParentGenos) <- c("P1", "P2", "P3", "P4") # Give names to the parents so that they
                                                # can be plotted easier

# The goal is to create three plots, one plot for each trait 
# (Oil, protein, SDWT) which show the distributions of each trait
# as a frequency histogram. There will be two histograms on each plot,
# one for each population that will be differentiated by the 
# fill of the bars. Each plot will also have the value of the populations'
# parents shown with a labeled arrow. 

## Section: Data wrangling
##################################################

# Pivot to a long format for easier plotting
AllPheno_long <- AllPheno %>%
  select(-Environment) %>%
  mutate(Oil     = Oil*10, 
         Protein = Protein*10) %>% # Convert to g/kg
  pivot_longer(c(SDWT, Oil, Protein)) %>%
  mutate(Population = factor(Population))

# Get the dataset for the parents 
AllPheno_long_parents <- AllPheno_long %>%
  filter(Genotype %in% ParentGenos)

AllPheno_long_parents$GenoShortName <- names(ParentGenos)[match(AllPheno_long_parents$Genotype, ParentGenos)]

# Split each dataframe by phenotype
AllPheno_long         <- split(AllPheno_long, AllPheno_long$name)
AllPheno_long_parents <- split(AllPheno_long_parents, AllPheno_long_parents$name)


## Section: Plotting functions
##################################################

# A dataframe to hold the units which a given trait is measured in
TraitUnits <- data.frame(Trait = c("Oil", "Protein", "SDWT"),
                         Units = c("Oil g/Kg^-1", "Protein g/Kg^-1", "g per 100 seeds"))

# The base plotting function. Takes as an argument a dataframe 
# containing measurements for a single phenotype (like AllPheno_long$oil)
Plot_initial <- function(TraitData = AllPheno_long$Oil, ParentData = AllPheno_long_parents$Oil, TraitName = "Oil"){
  
    if(TraitName == "Oil"){
      XLabel <- bquote("Oil g/"~Kg^-1)
    }else if(TraitName == "Protein"){
      XLabel <- bquote("Protein g/"~Kg^-1)
    }else{
      XLabel <- "grams per 100 seeds"
    }
  
  

  # Get the initial plot
  Plot.init <-   ggplot(TraitData, aes(x = value, fill = Population)) + 
    geom_histogram(position = position_dodge2(preserve = "single"), color = "black", alpha = 0.6, bins = 10) + 
    scale_fill_manual(values = c("#bdbdbd", "#636363")) + 
    theme_bw() + 
    ylab("Count") + 
    xlab(XLabel) + 
    theme_ipsum() + 
    theme(axis.text.x  = element_text(face = "bold", size = 20),
          axis.text.y  = element_text(size = 20),
          axis.title.y = element_text(face = "bold", size = 25),
          axis.title.x = element_text(face = "bold", size = 25),
          panel.border = element_rect(colour = "black", fill = NA),
          legend.title = element_text(face = "bold", size = 20),
          legend.text  = element_text(face = "bold", size = 20)) 
  
  # Data from the initial plot (I want the bin heights)
  PlotData <- ggplot_build(Plot.init)$data[[1]]
  
  # Using the bin counts from the plot, find the y-value where the labels for each
  # check/parent genotype shoud start
  ParentData$yval <- NA
  for(i in 1:nrow(ParentData)){
    ParentData$yval[[i]] <- PlotData$count[[max(which(PlotData$xmin < ParentData$value[[i]]))]]
  }
  
  # Add labels w/arrows for the parents/checks using this new data
  Plot.init + 
    ggrepel::geom_label_repel(data = ParentData,
                              aes(x = value, y = yval, label = GenoShortName),
                              nudge_y            = max(PlotData$count)/6,
                              arrow              = arrow(length = unit(0.015, "npc")),
                              min.segment.length = 1,
                              size               = 7,
                              show.legend        = FALSE,
                              inherit.aes        = FALSE)
}

# Make a plot for each trait
AllPlots        <- vector("list", length = length(AllPheno_long))
names(AllPlots) <- names(AllPheno_long)
for(i in 1:length(AllPlots)){
  CurrentTrait <- names(AllPlots)[[i]]
  AllPlots[[CurrentTrait]] <- Plot_initial(TraitData  = AllPheno_long[[CurrentTrait]], 
                                           ParentData = AllPheno_long_parents[[CurrentTrait]], 
                                           TraitName  = CurrentTrait)
}


# Lay out all three plots into a single image
FinalPlot <- plot_grid(plotlist   = AllPlots, 
                       labels     = c("A", "B", "C"),
                       align      = "v",
                       ncol       = 1,
                       label_size = 20)


# And save to an image
save_plot(paste0(here(), "/Plots/PhenoPlots.svg"),
          plot = FinalPlot,
          base_height = 13,
          base_width = 9)
