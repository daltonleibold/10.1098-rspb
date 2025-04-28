# source.R

rm(list = ls())
library(here)
library(tidyverse)
# load data
data <- read_csv(here("data", "data.csv"))
# source models from model comparisons
# - only run on initialization
# source(here("r", "model_comparisons.R))
# source summary statistics of raw data
source(here("r", "summary_stats.R"))
# source hypothesis testing output
source(here("r", "hypothesis_tests.R"))
# source plots
source(here("r", "plot_making.R"))

# unload everything we don't need
rm(list = c("brm.eggox.ld",
               "brm.eggox.ld.sup",
               "brm.eggox.lg",
               "brm.eggox.lg.sup",
               "brm.hmorph.ld",
               "brm.hmorph.lg",
               "fig.eggox",
               "fig.hmass",
               "fig.hsvl",
               "plot.data",
               "plot.eggox.ld",
               "plot.eggox.lg",
               "plot.hmass.ld",
               "plot.hmass.lg",
               "plot.hsvl.ld",
               "plot.hsvl.lg",
               "raw.data",
               "raw.eggox",
               "raw.hmass",
               "raw.hsvl"))
