# plots.R
# script for making plots

# load libraries
pacman::p_load(tidyverse, # for data wrangling
               brms, # for bayesian regression models 
               rstan, # for bayesian regression models
               flextable, # for making flextables
               svglite, # for making svg files
               ggpubr, # for arranging separate ggplots
               ggdist, # for shadeable density slabs
               gghalves, # for half-half geoms
               ggpp, # for position_dodge2nudge
               cowplot, # for publication-ready themes
               colorspace, # for lightening color palettes
               gridExtra, # for grid.arrange
               here) # for finding your files
# make a colorblind palette for plots

# set working directory
setwd(here())

# source custom functions
source(here("r", "func.r"))

# code for making custom plots
custom_theme <- theme_classic(base_size = 24) +
  #adjust axis title position
  theme(axis.title.y=element_text(vjust=1.5), 
        axis.title.x=element_text(vjust=0.2)) + 
  #adjust plot margins and line element size
  theme(plot.margin = unit(c(.3,.3,.6,.6), "cm"),
        line = element_line(size = 1.25)) + 
  #draw x and y axes
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) + 
  #put margins around axis labels so that nothing overlaps
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) + 
  # move tickmarks inside the axes and paint black
  theme(axis.ticks.length =unit(-0.3, "cm")) + 
  #spread out facets
  theme(panel.spacing = unit(2, units = "lines")) + 
  #make tick marks black
  theme(axis.ticks = element_line(color = "black")) + 
  #remove border from facet labels
  theme(strip.background = element_blank()) 

# make a colorblind palette for plots
cbPalette <- c("#ed7d31", "#f8cbad", 
               "#4472c4", "#b4c7e7")
cbPalette2 <- c("#b4c7e7", "#4472c4", 
                "#f8cbad", "#ed7d31")

# read in the raw data for making raw figures
data <- read_csv(here("data", "data.csv"))

# source models for making model figures
brm.eggox.ld <- readRDS(here("output", "models", "eggox_ld.rds"))
brm.eggox.lg <- readRDS(here("output", "models", "eggox_lg.rds"))
brm.hmorph.ld <- readRDS(here("output", "models", "hmorph_ld.rds"))
brm.hmorph.lg <- readRDS(here("output", "models", "hmorph_lg.rds"))

# extract posterior distributions from each model
# we'll use the models that included all independent effects of cort, temp, and egg mass
# rationale is that these help demonstrate the non-significant tendencies that aren't captured in the models of best fit
# for each of the models, we'll extract the posteriors of the main effects terms, excluding egg mass. we'll use these to calculate the predicted posterior distributions of each treatment level for the respective trait. for the multivariate model, we'll extract mass and svl separately. 

.plot.eggox.ld <- brm.eggox.ld %>%
  # extract the posterior distributions
  posterior_samples() %>%
  # select only the main effects terms
  select(grep("^b_", colnames(.))) %>%
  # add in trait and species columns
  # calculate predicted posterior distributions for each treatment level
  reframe(trait = "mortality",
          species = "delicata",
          A_23 = b_Intercept,
          B_23 = b_Intercept + b_cortB,
          A_28 = b_Intercept + b_temp28,
          B_28 = b_Intercept + b_cortB + b_temp28) %>%
  # pivot to long form
  pivot_longer(cols = A_23:B_28,
               names_to = "trt",
               values_to = "value") 

.plot.eggox.lg <- brm.eggox.lg %>%
  # extract the posterior distributions
  posterior_samples() %>%
  # select only the main effects terms
  select(grep("^b_", colnames(.))) %>%
  # add in trait and species columns
  # calculate predicted posterior distributions for each treatment level
  reframe(trait = "mortality",
          species = "guichenoti",
          A_23 = b_Intercept,
          B_23 = b_Intercept + b_cortB,
          A_28 = b_Intercept + b_temp28,
          B_28 = b_Intercept + b_cortB + b_temp28) %>%
  # pivot to long form
  pivot_longer(cols = A_23:B_28,
               names_to = "trt",
               values_to = "value")

.plot.hmass.ld <- brm.hmorph.ld %>%
  # extract the posterior distributions
  posterior_samples() %>%
  # select only the main effects terms
  # and only for hatchling mass
  select(grep("^b_hatchmass", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "hatchmass_")) %>%
  # add in trait and species columns
  # calculate predicted posterior distributions for each treatment level
  reframe(trait = "mass",
          species = "delicata",
          A_23 = b_Intercept,
          B_23 = b_Intercept + b_cortB,
          A_28 = b_Intercept + b_temp28,
          B_28 = b_Intercept + b_cortB + b_temp28) %>%
  # pivot to long form
  pivot_longer(cols = A_23:B_28,
               names_to = "trt",
               values_to = "value")

.plot.hmass.lg <- brm.hmorph.lg %>%
  # extract the posterior distributions
  posterior_samples() %>%
  # select only the main effects terms
  # and only for hatchling mass
  select(grep("^b_hatchmass", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "hatchmass_")) %>%
  # add in trait and species columns
  # calculate predicted posterior distributions for each treatment level
  reframe(trait = "mass",
          species = "guichenoti",
          A_23 = b_Intercept,
          B_23 = b_Intercept + b_cortB,
          A_28 = b_Intercept + b_temp28,
          B_28 = b_Intercept + b_cortB + b_temp28) %>%
  # pivot to long form
  pivot_longer(cols = A_23:B_28,
               names_to = "trt",
               values_to = "value")

.plot.hsvl.ld <- brm.hmorph.ld %>%
  # extract the posterior distributions
  posterior_samples() %>%
  # select only the main effects terms
  # and only for hatchling svl
  select(grep("^b_hatchsvl", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "hatchsvl_")) %>%
  # add in trait and species columns
  # calculate predicted posterior distributions for each treatment level
  reframe(trait = "svl",
          species = "delicata",
          A_23 = b_Intercept,
          B_23 = b_Intercept + b_cortB,
          A_28 = b_Intercept + b_temp28,
          B_28 = b_Intercept + b_cortB + b_temp28) %>%
  # pivot to long form
  pivot_longer(cols = A_23:B_28,
               names_to = "trt",
               values_to = "value")

.plot.hsvl.lg <- brm.hmorph.lg %>%
  # extract the posterior distributions
  posterior_samples() %>%
  # select only the main effects terms
  # and only for hatchling svl
  select(grep("^b_hatchsvl", colnames(.))) %>%
  rename_with(.cols = everything(),
              .fn = ~str_remove(.x, "hatchsvl_")) %>%
  # add in trait and species columns
  # calculate predicted posterior distributions for each treatment level
  reframe(trait = "svl",
          species = "guichenoti",
          A_23 = b_Intercept,
          B_23 = b_Intercept + b_cortB,
          A_28 = b_Intercept + b_temp28,
          B_28 = b_Intercept + b_cortB + b_temp28) %>%
  # pivot to long form
  pivot_longer(cols = A_23:B_28,
               names_to = "trt",
               values_to = "value")

# bind all the extracted posteriors together
plot.data <- bind_rows(.plot.eggox.ld,
                       .plot.eggox.lg,
                       .plot.hmass.ld,
                       .plot.hmass.lg,
                       .plot.hsvl.ld,
                       .plot.hsvl.lg) %>%
  # make a tracer column for pivoting
  group_by(trait, species, trt) %>%
  mutate(tracer = row_number()) %>%
  ungroup() %>%
  # pivot to wide form
  pivot_wider(names_from = "trait",
              values_from = "value") %>%
  # no longer need the tracer column
  select(-tracer) %>%
  # separate the trt column into its components
  separate(col = trt,
           into = c("cort", "temp"),
           sep = "_") %>%
  # recode everything to better names
  mutate(species = re.code(species, c("delicata" = "L. delicata",
                                      "guichenoti" = "L. guichenoti")),
         cort = re.code(cort, c("A" = "Control",
                                "B" = "CORT")),
         temp = re.code(temp, c("23" = "23°C",
                                "28" = "28°C"))) %>%
  # reunite the treatment column with the fixed names
  unite("trt", c(cort, temp), sep = ",") %>%
  # coerce species and treatment to factors
  mutate(species = factor(species, levels = c("L. delicata", "L. guichenoti")),
         trt = factor(trt, levels = c("CORT,28°C", "Control,28°C",
                                      "CORT,23°C", "Control,23°C")))



# plot embryonic mortality
fig.eggox <- plot.data %>%
  ggplot(aes(y = (exp(mortality)/(1+exp(mortality)))*100, 
             x = trt,
             fill = trt)) +   
  theme_half_open() +  
  scale_colour_manual(values = cbPalette, 
                      aesthetics = c("colour", "fill")) +  
  guides(fill_ramp = "none",
         color = guide_legend(override.aes = list(size = 5))) +
  stat_slab(side = "right", 
            scale = 0.4,
            show.legend = F,
            position = position_dodge(width = .8),
            aes(fill_ramp = stat(level)),.width = c(.95, .5, 1)) +
  scale_fill_ramp_discrete(from="gray95",
                           aesthetics = "fill_ramp") +
  geom_boxplot(width = .05,
               alpha = .5,
               outlier.alpha=0,
               position = position_dodge(width  = .8),
               show.legend = FALSE) +
  coord_flip() +
  facet_grid(cols = vars(species), 
             scales = "fixed") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  ylab("Egg Mortality (%)") +
  xlab(NULL) +
  custom_theme

# plot hatchling mass
fig.hmass <- plot.data %>%
  ggplot(aes(y = mass, 
             x = trt,
             fill = trt)) +   
  theme_half_open() +  
  scale_colour_manual(values = cbPalette, 
                      aesthetics = c("colour", "fill")) +  
  guides(fill_ramp = "none",
         color = guide_legend(override.aes = list(size = 5))) +
  stat_slab(side = "right", 
            scale = 0.4,
            show.legend = F,
            position = position_dodge(width = .8),
            aes(fill_ramp = stat(level)),.width = c(.95, .5, 1)) +
  scale_fill_ramp_discrete(from="gray95",
                           aesthetics = "fill_ramp") +
  geom_boxplot(width = .05,
               alpha = .5,
               outlier.alpha=0,
               position = position_dodge(width  = .8),
               show.legend = FALSE) +
  coord_flip() +
  facet_grid(cols = vars(species), 
             scales = "fixed") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  ylab("Hatchling Mass (g)") +
  xlab(NULL) +
  custom_theme

# plot hatchling svl
fig.hsvl <- plot.data %>%
  ggplot(aes(y = svl, 
             x = trt,
             fill = trt)) +   
  theme_half_open() +  
  scale_colour_manual(values = cbPalette, 
                      aesthetics = c("colour", "fill")) +  
  guides(fill_ramp = "none",
         color = guide_legend(override.aes = list(size = 5))) +
  stat_slab(side = "right", 
            scale = 0.4,
            show.legend = F,
            position = position_dodge(width = .8),
            aes(fill_ramp = stat(level)),.width = c(.95, .5, 1)) +
  scale_fill_ramp_discrete(from="gray95",
                           aesthetics = "fill_ramp") +
  geom_boxplot(width = .05,
               alpha = .5,
               outlier.alpha=0,
               position = position_dodge(width  = .8),
               show.legend = FALSE) +
  coord_flip() +
  facet_grid(cols = vars(species), 
             scales = "fixed") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  ylab("Hatchling SVL (mm)") +
  xlab(NULL) +
  custom_theme

# prep the raw data for plotting
raw.data <- data %>%
  mutate(species = re.code(species, c("delicata" = "L. delicata",
                                      "guichenoti" = "L. guichenoti")),
         cort = re.code(cort, c("A" = "Control",
                                "B" = "CORT")),
         temp = re.code(temp, c("23" = "23°C",
                                "28" = "28°C"))) %>%
  unite("trt", c(cort, temp), sep = ",\n") %>%
  mutate(species = factor(species,
                          levels = c("L. delicata", "L. guichenoti")),
         trt = factor(trt, levels = c("Control,\n23°C", "CORT,\n23°C",
                                      "Control,\n28°C", "CORT,\n28°C")))

# plot hatchling mortality
raw.eggox <- raw.data %>%
  filter(!is.na(egg_death)) %>%
  ggplot(aes(x = trt, fill = interaction(egg_death, trt))) +
  geom_bar(position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(cols = vars(species)) +
  scale_fill_manual(values = c("black", "#b4c7e7", 
                               "black", "#4472c4",
                               "black", "#f8cbad", 
                               "black", "#ed7d31")) +
  guides(fill = "none") +
  xlab(NULL) +
  ylab("Hatching Success") +
  custom_theme

# plot hatchling mass
raw.hmass <- raw.data %>%
  ggplot(aes(x = trt, y = hatch_mass, col = trt)) +
  scale_color_manual(values = cbPalette2) +
  stat_summary() +
  guides(color = FALSE) +
  facet_wrap(~species, scales = "free_y") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  xlab(NULL) +
  ylab("Hatchling Mass (g)") +
  custom_theme

# plot hatchling svl
raw.hsvl <- raw.data %>%
  ggplot(aes(x = trt, y = hatch_svl, col = trt)) +
  scale_color_manual(values = cbPalette2) +
  stat_summary() +
  guides(color = FALSE) +
  facet_wrap(~species, scales = "free_y") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  xlab(NULL) +
  ylab("Hatchling SVL (mm)") +
  custom_theme


ggsave(plot = fig.eggox,
       filename = "fig_2_a_i - model egg mortality.png", 
       path = here("output", "figures", "fig2 - data plots"),
       width = 3600,
       height = 2400,
       units = "px")
ggsave(plot = raw.eggox,
       filename = "fig_2_a_ii - raw egg mortality.png", 
       path = here("output", "figures", "fig2 - data plots"),
       width = 3600,
       height = 2400,
       units = "px")
ggsave(plot = fig.hmass,
       filename = "fig_2_b_i - model hatchling mass.png", 
       path = here("output", "figures", "fig2 - data plots"),
       width = 3600,
       height = 2400,
       units = "px")
ggsave(plot = raw.hmass,
       filename = "fig_2_b_ii - raw hatchling mass.png", 
       path = here("output", "figures", "fig2 - data plots"),
       width = 3600,
       height = 2400,
       units = "px")
ggsave(plot = fig.hsvl,
       filename = "fig_2_c_i - model hatchling svl.png", 
       path = here("output", "figures", "fig2 - data plots"),
       width = 3600,
       height = 2400,
       units = "px")
ggsave(plot = raw.hsvl,
       filename = "fig_2_c_ii - raw hatchling svl.png", 
       path = here("output", "figures", "fig2 - data plots"),
       width = 3600,
       height = 2400,
       units = "px")

rm(raw.data)

