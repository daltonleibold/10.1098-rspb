# model_comparisons.R
# script for constructing and comparing models

#----
# setup
#----
# clear workspace
base::rm(list = ls())
# set relative working directory
pacman::p_load(here)
base::setwd(here::here())
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
               gridExtra) # for grid.arrange
# source custom functions
base::source(here("r", "func.R"))
# make a colorblind palette for plots
cbPalette <- c("#b4c7e7", "#f8cbad",
               "#4472c4", "#ed7d31")
#----
# data prep
#----
# read in data
data <- read_csv(here("data", "data.csv"))
# prepare data for modelling
.data.mod <- data %>%
  # center the morphology variables about the species-level means
  group_by(species) %>%
  mutate(across(.cols = c(hatch_mass:hatch_tl, 
                          egg_mass:egg_width),
                .fns = ~scale(.x, center = T, scale = F))) %>%
  ungroup() %>%
  # parse columns to correct formats
  mutate(across(.cols = c(species, egg_id, clutch),
                .fns = factor),
         cort = factor(cort, levels = c("A", "B")),
         temp = factor(temp, levels = c("23", "28")),
         egg_death = factor(egg_death, levels = c("hatched", "died")),
         egg_ox = as.logical(egg_ox)) %>%
  # remove any animals we did not measure mass on
  dplyr::filter(!is.na(egg_mass)) 
# take a look
glimpse(.data.mod)

# dependent vars
## egg mortality
## hatchling mass
## hatchling SVL

# independent vars
## cort
## temp
## egg mass

# grouping vars
## species

# random effects
## clutch

##########################
# Model Comparisons
##########################
#----
# embryonic mortality
#----
## possible model structures
## each permutation of cort, temp, and egg mass
## with clutch id included as a random effect
# additive cort, temp, egg_mass
.m1 <- bf(egg_ox ~ cort + temp + egg_mass + (1|clutch)) + bernoulli(link = "logit")
# interaction between cort and temp with independent effect of egg mass
.m2 <- bf(egg_ox ~ cort * temp + egg_mass + (1|clutch)) + bernoulli(link = "logit")
# interaction between cort and egg mass with independent effect of temp
.m3 <- bf(egg_ox ~ cort * egg_mass + temp + (1|clutch)) + bernoulli(link = "logit")
# interaction between temp and egg mass with independent effect of cort
.m4 <- bf(egg_ox ~ temp * egg_mass + cort + (1|clutch)) + bernoulli(link = "logit")
# all two way interactions of cort, temp, and egg mass
.m5 <- bf(egg_ox ~ (cort + temp + egg_mass)^2 + (1|clutch)) + bernoulli(link = "logit")
# all three-way interaction of cort, temp, and egg mass
.m6 <- bf(egg_ox ~ cort * temp * egg_mass + (1|clutch)) + bernoulli(link = "logit")
# make a list of the possible models
.mods <- list(.m1, .m2, .m3,
              .m4, .m5, .m6)
#----
## delicata
#----
# make an empty list to hold your brmsfit objects
.brm <- list()
# for each formula structure, run a simple bayesian regression model and store it
## we model each species separately
## we can compare the distributions among species anyways
for(i in 1:length(.mods)){
  .brm[[i]] <- brm(.mods[[i]], data = .data.mod[.data.mod$species == "delicata",])
}
# fit the loo criterion to each model
.fit <- lapply(.brm, function(x){add_criterion(x, criterion = "loo")})
# compare the models
.loo_comp <- loo_compare(.fit[[1]], .fit[[2]], 
                         .fit[[3]], .fit[[4]], 
                         .fit[[5]], .fit[[6]],
                         model_names = c("Egg Mortality ~ CORT + Temp + Egg Mass + (1|Clutch)",
                                         "Egg Mortality ~ CORT * Temp + Egg Mass + (1|Clutch)",
                                         "Egg Mortality ~ CORT * Egg Mass + Temp + (1|Clutch)",
                                         "Egg Mortality ~ CORT + Temp * Egg Mass + (1|Clutch)",
                                         "Egg Mortality ~ (CORT + Temp + Egg Mass)^2 + (1|Clutch)",
                                         "Egg Mortality ~ CORT * Temp * Egg Mass + (1|Clutch)"))
# take a look
# looks like the best model is model 1: egg ox ~ cort + temp + egg mass
tbl.s.4.a.i.ld <- .loo_comp %>%
  as_tibble() %>%
  mutate(Model = row.names(.loo_comp),
         Species = "L. delicata") %>%
  mutate(elpd_diff = round(elpd_diff, digits = 3),
         se_diff = round(se_diff, digits = 3)) %>%
  mutate(`Δ ELPD ± S.E.` = paste(elpd_diff, se_diff, sep = " ± ")) %>%
  select(Species, Model, `Δ ELPD ± S.E.`) %>%
  flextable() %>%
  width(j = 1, width = 1.5) %>%
  width(j = 2, width = 5) %>%
  width(j = 3, width = 2) %>%
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top") %>%
  bg(bg = "white", part = "all")
tbl.s.4.a.i.ld

# save the flextables 
save_as_image(tbl.s.4.a.i.ld, here("output", "tables", "tbl_s4", "tbl_s4_a_i_ld.png"))

# best formula:
.bf.eggox.ld <- bf(egg_ox ~ cort + temp + egg_mass + (1|clutch)) + bernoulli(link = "logit")

#----
## guichenoti
#----
# make an empty list to hold your brmsfit objects
.brm <- list()
# for each formula structure, run a simple bayesian regression model and store it
for(i in 1:length(.mods)){
  .brm[[i]] <- brm(.mods[[i]], data = .data.mod[.data.mod$species == "guichenoti",])
}
# fit the loo criterion to each model
.fit <- lapply(.brm, function(x){add_criterion(x, criterion = "loo")})
# compare the models
.loo_comp <- loo_compare(.fit[[1]], .fit[[2]], .fit[[3]],
                         .fit[[4]], .fit[[5]], .fit[[6]],
                         model_names = c("Egg Mortality ~ CORT + Temp + Egg Mass + (1|Clutch)",
                                         "Egg Mortality ~ CORT * Temp + Egg Mass + (1|Clutch)",
                                         "Egg Mortality ~ CORT * Egg Mass + Temp + (1|Clutch)",
                                         "Egg Mortality ~ CORT + Temp * Egg Mass + (1|Clutch)",
                                         "Egg Mortality ~ (CORT + Temp + Egg Mass)^2 + (1|Clutch)",
                                         "Egg Mortality ~ CORT * Temp * Egg Mass + (1|Clutch)"))
# take a look
# looks like all of the models are identical, presumably driven by the high variance caused by the overall low mortality in guichenoti.
# we'll compare all of the models including clutch as a random effect and go from there
tbl.s.4.a.i.lg <- .loo_comp %>%
  as_tibble() %>%
  mutate(Model = row.names(.loo_comp),
         Species = "L. guichenoti") %>%
  mutate(elpd_diff = round(elpd_diff, digits = 3),
         se_diff = round(se_diff, digits = 3)) %>%
  mutate(`Δ ELPD ± S.E.` = paste(elpd_diff, se_diff, sep = " ± ")) %>%
  select(Species, Model, `Δ ELPD ± S.E.`) %>%
  flextable() %>%
  width(j = 1, width = 1.5) %>%
  width(j = 2, width = 5) %>%
  width(j = 3, width = 2) %>%
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top") %>%
  bg(bg = "white", part = "all")
tbl.s.4.a.i.lg


# the model including all three way interactions is *technically* the best (lowest dELPD), but the error is extremely high in all other models, and they are statistically indistinguishable from one another. 
# with how few observations of mortality we have, i suspect a three-way interaction is overfitting. 
# since no models including interaction terms are any more or less informative than the additive null model, we will default to using the additive null model.

# best formula:
.bf.eggox.lg <- bf(egg_ox ~ cort + temp + egg_mass + (1|clutch)) + bernoulli(link = "logit")

# save the flextables
save_as_image(tbl.s.4.a.i.lg, here("output", "tables", "tbl_s4", "tbl_s4_a_i_lg.png"))
#----
# hatchling mass
#----
## possible model structures
# additive cort, temp, egg_mass
.m1 <- bf(hatch_mass ~ cort + temp + egg_mass + (1|clutch)) + gaussian()
# interaction between cort and temp with independent effect of egg mass
.m2 <- bf(hatch_mass ~ cort * temp + egg_mass + (1|clutch)) + gaussian()
# interaction between cort and egg mass with independent effect of temp
.m3 <- bf(hatch_mass ~ cort * egg_mass + temp + (1|clutch)) + gaussian()
# interaction between temp and egg mass with independent effect of cort
.m4 <- bf(hatch_mass ~ temp * egg_mass + cort + (1|clutch)) + gaussian()
# all two way interactions of cort, temp, and egg mass
.m5 <- bf(hatch_mass ~ (cort + temp + egg_mass)^2 + (1|clutch)) + gaussian()
# all three-way interaction of cort, temp, and egg mass
.m6 <- bf(hatch_mass ~ cort * temp * egg_mass + (1|clutch)) + gaussian()
# make a list of the possible models
.mods <- list(.m1, .m2, .m3,
              .m4, .m5, .m6)
#----
## delicata
#----
# make an empty list to hold your brmsfit objects
.brm <- list()
# for each formula structure, run a simple bayesian regression model and store it
for(i in 1:length(.mods)){
  .brm[[i]] <- brm(.mods[[i]], data = .data.mod[.data.mod$species == "delicata",])
}
# fit the loo criterion to each model
.fit <- lapply(.brm, function(x){add_criterion(x, criterion = "loo")})
# compare the models
.loo_comp <- loo_compare(.fit[[1]], .fit[[2]], .fit[[3]],
                         .fit[[4]], .fit[[5]], .fit[[6]],
                         model_names = c("Hatchling Mass ~ CORT + Temp + Egg Mass + (1|Clutch)",
                                         "Hatchling Mass ~ CORT * Temp + Egg Mass + (1|Clutch)",
                                         "Hatchling Mass ~ CORT * Egg Mass + Temp + (1|Clutch)",
                                         "Hatchling Mass ~ CORT + Temp * Egg Mass + (1|Clutch)",
                                         "Hatchling Mass ~ (CORT + Temp + Egg Mass)^2 + (1|Clutch)",
                                         "Hatchling Mass ~ CORT * Temp * Egg Mass + (1|Clutch)"))
# take a look
tbl.s.4.b.i.ld <- .loo_comp %>%
  as_tibble() %>%
  mutate(Model = row.names(.loo_comp),
         Species = "L. delicata") %>%
  mutate(elpd_diff = round(elpd_diff, digits = 3),
         se_diff = round(se_diff, digits = 3)) %>%
  mutate(`Δ ELPD ± S.E.` = paste(elpd_diff, se_diff, sep = " ± ")) %>%
  select(Species, Model, `Δ ELPD ± S.E.`) %>%
  flextable() %>%
  width(j = 1, width = 1.5) %>%
  width(j = 2, width = 5) %>%
  width(j = 3, width = 2) %>%
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top") %>%
  bg(bg = "white", part = "all")
tbl.s.4.b.i.ld
# save the flextable
save_as_image(tbl.s.4.b.i.ld, here("output", "tables", "tbl_s4", "tbl_s4_b_i_ld.png"))
# best formula:
# once again, the models are indistinguishable, so we opt
# to use the additive null model
.bf.mass.ld <- bf(hatch_mass ~ cort + temp + egg_mass + (1|clutch)) + gaussian()

#----
## guichenoti
#----
# make an empty list to hold your brmsfit objects
.brm <- list()
# for each formula structure, run a simple bayesian regression model and store it
for(i in 1:length(.mods)){
  .brm[[i]] <- brm(.mods[[i]], data = .data.mod[.data.mod$species == "guichenoti",])
}
# fit the loo criterion to each model
.fit <- lapply(.brm, function(x){add_criterion(x, criterion = "loo")})
# compare the models
.loo_comp <- loo_compare(.fit[[1]], .fit[[2]], .fit[[3]],
                         .fit[[4]], .fit[[5]], .fit[[6]],
                         model_names = c("Hatchling Mass ~ CORT + Temp + Egg Mass + (1|Clutch)",
                                         "Hatchling Mass ~ CORT * Temp + Egg Mass + (1|Clutch)",
                                         "Hatchling Mass ~ CORT * Egg Mass + Temp + (1|Clutch)",
                                         "Hatchling Mass ~ CORT + Temp * Egg Mass + (1|Clutch)",
                                         "Hatchling Mass ~ (CORT + Temp + Egg Mass)^2 + (1|Clutch)",
                                         "Hatchling Mass ~ CORT * Temp * Egg Mass + (1|Clutch)"))
# take a look
tbl.s.4.b.i.lg <- .loo_comp %>%
  as_tibble() %>%
  mutate(Model = row.names(.loo_comp),
         Species = "L. guichenoti") %>%
  mutate(elpd_diff = round(elpd_diff, digits = 3),
         se_diff = round(se_diff, digits = 3)) %>%
  mutate(`Δ ELPD ± S.E.` = paste(elpd_diff, se_diff, sep = " ± ")) %>%
  select(Species, Model, `Δ ELPD ± S.E.`) %>%
  flextable() %>%
  width(j = 1, width = 1.5) %>%
  width(j = 2, width = 5) %>%
  width(j = 3, width = 2) %>%
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top") %>%
  bg(bg = "white", part = "all")
tbl.s.4.b.i.lg
# save the flextable
save_as_image(tbl.s.4.b.i.lg, here("output", "tables", "tbl_s4", "tbl_s4_b_i_lg.png"))
# best formula:
.bf.mass.lg <- bf(hatch_mass ~ cort + temp + egg_mass + (1|clutch)) + gaussian()

#----
# hatchling SVL
#----
## possible model structures
# additive cort, temp, egg_mass
.m1 <- bf(hatch_svl ~ cort + temp + egg_mass + (1|clutch)) + gaussian()
# interaction between cort and temp with independent effect of egg mass
.m2 <- bf(hatch_svl ~ cort * temp + egg_mass + (1|clutch)) + gaussian()
# interaction between cort and egg mass with independent effect of temp
.m3 <- bf(hatch_svl ~ cort * egg_mass + temp + (1|clutch)) + gaussian()
# interaction between temp and egg mass with independent effect of cort
.m4 <- bf(hatch_svl ~ temp * egg_mass + cort + (1|clutch)) + gaussian()
# all two way interactions of cort, temp, and egg mass
.m5 <- bf(hatch_svl ~ (cort + temp + egg_mass)^2 + (1|clutch)) + gaussian()
# all three-way interaction of cort, temp, and egg mass
.m6 <- bf(hatch_svl ~ cort * temp * egg_mass + (1|clutch)) + gaussian()
# make a list of the possible models
.mods <- list(.m1, .m2, .m3,
              .m4, .m5, .m6)
#----
## delicata
#----
# make an empty list to hold your brmsfit objects
.brm <- list()
# for each formula structure, run a simple bayesian regression model and store it
for(i in 1:length(.mods)){
  .brm[[i]] <- brm(.mods[[i]], data = .data.mod[.data.mod$species == "delicata",])
}
# fit the loo criterion to each model
.fit <- lapply(.brm, function(x){add_criterion(x, criterion = "loo")})
# compare the models
.loo_comp <- loo_compare(.fit[[1]], .fit[[2]], .fit[[3]],
                         .fit[[4]], .fit[[5]], .fit[[6]],
                         model_names = c("Hatchling SVL ~ CORT + Temp + Egg Mass + (1|Clutch)",
                                         "Hatchling SVL ~ CORT * Temp + Egg Mass + (1|Clutch)",
                                         "Hatchling SVL ~ CORT * Egg Mass + Temp + (1|Clutch)",
                                         "Hatchling SVL ~ CORT + Temp * Egg Mass + (1|Clutch)",
                                         "Hatchling SVL ~ (CORT + Temp + Egg Mass)^2 + (1|Clutch)",
                                         "Hatchling SVL ~ CORT * Temp * Egg Mass + (1|Clutch)"))
# take a look
tbl.s.4.c.i.ld <- .loo_comp %>%
  as_tibble() %>%
  mutate(Model = row.names(.loo_comp),
         Species = "L. delicata") %>%
  mutate(elpd_diff = round(elpd_diff, digits = 3),
         se_diff = round(se_diff, digits = 3)) %>%
  mutate(`Δ ELPD ± S.E.` = paste(elpd_diff, se_diff, sep = " ± ")) %>%
  select(Species, Model, `Δ ELPD ± S.E.`) %>%
  flextable() %>%
  width(j = 1, width = 1.5) %>%
  width(j = 2, width = 5) %>%
  width(j = 3, width = 2) %>%
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top") %>%
  bg(bg = "white", part = "all")
tbl.s.4.c.i.ld
# save the flextable
save_as_image(tbl.s.4.c.i.ld, here("output", "tables", "tbl_s4", "tbl_s4_c_i_ld.png"))
# best formula:
.bf.mass.ld <- bf(hatch_svl ~ cort + temp + egg_mass + (1|clutch)) + gaussian()

#----
## guichenoti
#----
# make an empty list to hold your brmsfit objects
.brm <- list()
# for each formula structure, run a simple bayesian regression model and store it
for(i in 1:length(.mods)){
  .brm[[i]] <- brm(.mods[[i]], data = .data.mod[.data.mod$species == "guichenoti",])
}
# fit the loo criterion to each model
.fit <- lapply(.brm, function(x){add_criterion(x, criterion = "loo")})
# compare the models
.loo_comp <- loo_compare(.fit[[1]], .fit[[2]], .fit[[3]],
                         .fit[[4]], .fit[[5]], .fit[[6]],
                         model_names = c("Hatchling SVL ~ CORT + Temp + Egg Mass + (1|Clutch)",
                                         "Hatchling SVL ~ CORT * Temp + Egg Mass + (1|Clutch)",
                                         "Hatchling SVL ~ CORT * Egg Mass + Temp + (1|Clutch)",
                                         "Hatchling SVL ~ CORT + Temp * Egg Mass + (1|Clutch)",
                                         "Hatchling SVL ~ (CORT + Temp + Egg Mass)^2 + (1|Clutch)",
                                         "Hatchling SVL ~ CORT * Temp * Egg Mass + (1|Clutch)"))
# take a look
tbl.s.4.c.i.lg <- .loo_comp %>%
  as_tibble() %>%
  mutate(Model = row.names(.loo_comp),
         Species = "L. guichenoti") %>%
  mutate(elpd_diff = round(elpd_diff, digits = 3),
         se_diff = round(se_diff, digits = 3)) %>%
  mutate(`Δ ELPD ± S.E.` = paste(elpd_diff, se_diff, sep = " ± ")) %>%
  select(Species, Model, `Δ ELPD ± S.E.`) %>%
  flextable() %>%
  width(j = 1, width = 1.5) %>%
  width(j = 2, width = 5) %>%
  width(j = 3, width = 2) %>%
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top") %>%
  bg(bg = "white", part = "all")
tbl.s.4.c.i.lg
# save the flextable
save_as_image(tbl.s.4.c.i.lg, here("output", "tables", "tbl_s4", "tbl_s4_c_i_lg.png"))
# best formula:
.bf.svl.lg <- bf(hatch_svl ~ cort + temp + egg_mass + (1|clutch)) + gaussian()


#########################
# running best models
#########################

# now that we have run model comparisons and found the best model structures for each response variable for each species, we can now run the actual models.
# first, let's put all of the best models here for easy reference
## the best model for embryonic mortality includes the independent effects of cort treatment, temp treatment, and starting egg mass, but no interaction terms
.bf.eggox <- bf(egg_ox ~ cort + temp + egg_mass + (1|clutch)) + bernoulli(link = "logit")
## the best formula for mass and svl include the independent effects of cort treatment, temperature treatment, and egg mass, but no interaction terms.
.bf.mass <- bf(hatch_mass ~ cort + temp + egg_mass + (1|clutch)) + gaussian()
.bf.svl <- bf(hatch_svl ~ cort + temp + egg_mass + (1|clutch)) + gaussian()
## in all models, egg clutch is included as a random intercept to account for sibling effects

# let's run the univariate model for embryonic mortality
## for delicata first
.brm.eggox.ld <- brm(formula = .bf.eggox,
                     data = .data.mod[.data.mod$species == "delicata",],
                     cores = 4, chains = 4,
                     warmup = 2500, iter = 5000,
                     file_refit = "on_change",
                     file = here("output", "models", "eggox_ld")
)
## for guichenoti next
.brm.eggox.lg <- brm(formula = .bf.eggox,
                     data = .data.mod[.data.mod$species == "guichenoti",],
                     cores = 4, chains = 4,
                     warmup = 2500, iter = 5000,
                     file_refit = "on_change",
                     file = here("output", "models", "eggox_lg")
)

# now let's run the multivariate model for hatchling morphology
## for delicata first
.brm.hmorph.ld <- brm(formula = .bf.mass + .bf.svl + set_rescor(rescor = T),
                      data = .data.mod[.data.mod$species == "delicata",],
                      cores = 4, chains = 4,
                      warmup = 2500, iter = 5000,
                      file_refit = "on_change",
                      file = here("output", "models", "hmorph_ld")
)
## for guichenoti next
.brm.hmorph.lg <- brm(formula = .bf.mass + .bf.svl + set_rescor(rescor = T),
                      data = .data.mod[.data.mod$species == "guichenoti",],
                      cores = 4, chains = 4,
                      warmup = 2500, iter = 5000,
                      file_refit = "on_change",
                      file = here("output", "models", "hmorph_lg")
)

