# This main script calls all subscripts

# set working directory and empty environment
setwd("~/Projects/SLR/deep_uncertainty/analysis/R-scripts")
rm(list = ls())
graphics.off()

# load libraries for plotting and set plotting preferences
source("set_theme.R")
theme_set(theme_ali(base_size=20))

# Compare recent (local) assessments
source("plot_recent_assessments.R")

# Analyse affect of different interpretation IPCC's likely' range
source("likely_range.R")

