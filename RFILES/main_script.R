# This main script calls all subscripts
# Questions, contact Alexander Bakker: alexander.bakker@rws.nl
# ================================================================================================

# set working directory and empty environment
#setwd("~/Projects/SLR/deep_uncertainty/analysis/R-scripts")
rm(list = ls())
graphics.off()

# load libraries for plotting
library(reshape2)
library(lhs)
library(rriskDistributions)
library(ggplot2)
library(grid)
library(gridExtra)

# Compare recent (local) assessments
source("plot_recent_assessments.R")

# Analyse affect of different interpretation IPCC's likely' range
source("likely_range.R")

# Analyse affect of different interpretation on West-Antarctic ice sheet contribution to sea-level
source("consensus_assessment.R")

# NOTE: warnings will appear when plotting. This is because some rows contain NA's and therefore 
#       these warnings can be ignored. For example, Perrette et al. 2013 does not project land 
#       contribution to GSL.
