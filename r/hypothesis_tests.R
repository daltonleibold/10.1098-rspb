# model summaries, posterior distribution extraction, and pairwise comparisons
# - setup
#----
pacman::p_load(tidyverse, brms, rstan, flextable, here)
setwd(here::here())
set_flextable_defaults(
  font.family = "Arial",
  font.size = 12,
  digits = 3
  )
.fx.ft <- function(x){
  x %>%
    flextable() %>%
    padding(padding = 10, part = "all") %>%
    width(width = 1.5) %>%
    width(width = 3,
          j = 2) %>%
    width(width = 3,
          j = 3) %>%
    merge_v(j = 1:2) %>%
    valign(j = 1:2, valign = "top") %>%
    align(align = "center", part = "all") %>%
    bg(bg = "white", part = "all")
}
source(here::here("r", "func.R"))
# exp(ln(odds)) / 1 + (exp(ln(odds)))
#----
# - source the models
.brm.eggox.ld <- readRDS(here("output", "models", "eggox_ld.rds"))
.brm.eggox.lg <- readRDS(here("output", "models", "eggox_lg.rds"))
.brm.hmorph.ld <- readRDS(here("output", "models", "hmorph_ld.rds"))
.brm.hmorph.lg <- readRDS(here("output", "models", "hmorph_lg.rds"))
#----
# - make the model summaries
#----
# -- egg mortality
#----
# | best model: egg mortality ~ cort + temp + egg mass + (1|clutch)
# --- lampropholis delicata
.tbl.2.a.i <- mod.sum(.brm.eggox.ld) %>%
  mutate(Trait = "Egg Mortality",
         Species = "L. delicata",
         Contrast = c("Intercept", 
                      "Stress\n[CORT - Control]",
                      "Temperature\n[28°C - 23°C]",
                      "Egg Mass")) %>%
  select(Trait, Species, Contrast, everything()) 
# --- lampropholis guichenoti
.tbl.2.a.ii <- mod.sum(.brm.eggox.lg) %>%
  mutate(Trait = "Egg Mortality",
         Species = "L. guichenoti",
         Contrast = c("Intercept", 
                      "Stress\n[CORT - Control]",
                      "Temperature\n[28°C - 23°C]",
                      "Egg Mass")) %>%
  select(Trait, Species, Contrast, everything())
#----
# -- hatchling mass
#----
# -| best model: hatchling mass ~ cort + temp + egg mass + (1|clutch)
# --- lampropholis delicata
.tbl.2.b.i <- mod.sum(.brm.hmorph.ld) %>%
  mutate(Contrast = row.names(.)) %>%
  filter(str_detect(Contrast, "^hatchmass") == T) %>%
  mutate(Trait = "Hatchling Mass",
         Species = "L. delicata",
         Contrast = c("Intercept", 
                      "Stress\n[CORT - Control]",
                      "Temperature\n[28°C - 23°C]",
                      "Egg Mass")) %>%
  select(Trait, Species, Contrast, everything()) 
# --- lampropholis guichenoti
.tbl.2.b.ii <- mod.sum(.brm.hmorph.lg) %>%
  mutate(Contrast = row.names(.)) %>%
  filter(str_detect(Contrast, "^hatchmass") == T) %>%
  mutate(Trait = "Hatchling Mass",
         Species = "L. guichenoti",
         Contrast = c("Intercept", 
                      "Stress\n[CORT - Control]",
                      "Temperature\n[28°C - 23°C]",
                      "Egg Mass")) %>%
  select(Trait, Species, Contrast, everything())
#----
# -- hatchling svl
#----
# -| best model: hatchling mass ~ cort + temp + egg mass + (1|clutch)
# --- lampropholis delicata
.tbl.2.c.i <- mod.sum(.brm.hmorph.ld) %>%
  mutate(Contrast = row.names(.)) %>%
  filter(str_detect(Contrast, "^hatchsvl") == T) %>%
  mutate(Trait = "Hatchling SVL",
         Species = "L. delicata",
         Contrast = c("Intercept", 
                      "Stress\n[CORT - Control]",
                      "Temperature\n[28°C - 23°C]",
                      "Egg Mass")) %>%
  select(Trait, Species, Contrast, everything()) 
# --- lampropholis guichenoti
.tbl.2.c.ii <- mod.sum(.brm.hmorph.lg) %>%
  mutate(Contrast = row.names(.)) %>%
  filter(str_detect(Contrast, "^hatchsvl") == T) %>%
  mutate(Trait = "Hatchling SVL",
         Species = "L. guichenoti",
         Contrast = c("Intercept", 
                      "Stress\n[CORT - Control]",
                      "Temperature\n[28°C - 23°C]",
                      "Egg Mass")) %>%
  select(Trait, Species, Contrast, everything())
#----
# - compiling all the main stats
#----
# -| table 2 - stats tests (main effects)
# --- merge the model summaries
.tbl.2.a <- bind_rows(.tbl.2.a.i, .tbl.2.a.ii)
.tbl.2.b <- bind_rows(.tbl.2.b.i, .tbl.2.b.ii)
.tbl.2.c <- bind_rows(.tbl.2.c.i, .tbl.2.c.ii)
.tbl.2 <- bind_rows(.tbl.2.a, .tbl.2.b)
# --- make flextables
.tbl.2.a.i <- .fx.ft(.tbl.2.a.i)
.tbl.2.a.ii <- .fx.ft(.tbl.2.a.ii)
.tbl.2.b.i <- .fx.ft(.tbl.2.b.i)
.tbl.2.b.ii <- .fx.ft(.tbl.2.b.ii)
.tbl.2.c.i <- .fx.ft(.tbl.2.c.i)
.tbl.2.c.ii <- .fx.ft(.tbl.2.c.ii)
.tbl.2.a <- .fx.ft(.tbl.2.a)
.tbl.2.b <- .fx.ft(.tbl.2.b)
.tbl.2.c <- .fx.ft(.tbl.2.c)
.tbl.2 <- .fx.ft(.tbl.2)
# --- save the flextable for table 2
save_as_image(.tbl.2, here("output", "tables", "tbl2 - stats tests.png"))
stats.tests <- .tbl.2
#----
# - extract the posterior distributions from each model
#----
# -- egg mortality
#----
# | best model: egg mortality ~ cort + temp + egg mass + (1|clutch)
# --- lampropholis delicata
.post.eggox.ld <- .brm.eggox.ld %>%
  # extract the posterior distributions
  posterior_samples() %>%
  # select only the main effects (start with "b_")
  select(grep("^b_", colnames(.))) %>%
  # calculate predicted posterior distributions for each treatment level
  reframe(Trait = "Egg Mortality", # add a trait column
          Species = "L. delicata", # add a species column
          `Control,23°C` = b_Intercept, # Control,23 is the default
          `Control,28°C` = b_Intercept + b_temp28, # Add effect of temp
          `CORT,23°C` = b_Intercept + b_cortB, # add effect of cort
          `CORT,28°C` = b_Intercept + b_cortB + b_temp28) # add effect of temp and cort
# --- lampropholis guichenoti
.post.eggox.lg <- .brm.eggox.lg %>%
  # extract the posterior distributions
  posterior_samples() %>%
  # select only the main effects (start with "b_)
  select(grep("^b_", colnames(.))) %>%
  # calculate the predicted posterior distributions for each treatment level
  reframe(Trait = "Egg Mortality", # add a trait column
          Species = "L. guichenoti", # add a species column
          `Control,23°C` = b_Intercept, # Control, 23 is the default
          `Control,28°C` = b_Intercept + b_temp28, # add effect of temp
          `CORT,23°C` = b_Intercept + b_cortB, # add effect of cort
          `CORT,28°C` = b_Intercept + b_cortB + b_temp28) # add effect of temp and cort
#----
# -- hatchling mass
#----
# -| best model: hatchling mass ~ cort + temp + egg mass + (1|clutch)
# --- lampropholis delicata
.post.hmass.ld <- .brm.hmorph.ld %>%
  # extract posterior distribuitions
  posterior_samples() %>%
  # select only the main effects (start with "b_")
  select(grep("^b_hatchmass", colnames(.))) %>%
  # calculate the predicted posterior distributions for each treatment level
  reframe(Trait = "Hatchling Mass", # add a trait column
          Species = "L. delicata", # add a species column
          `Control,23°C` = b_hatchmass_Intercept, # Control, 23 is the default
          `Control,28°C` = b_hatchmass_Intercept + b_hatchmass_temp28, # add effect of temp
          `CORT,23°C` = b_hatchmass_Intercept + b_hatchmass_cortB, # add effect of cort
          `CORT,28°C` = b_hatchmass_Intercept + b_hatchmass_cortB + b_hatchmass_temp28) # add effect of temp and cort
# --- lampropholis guichenoti
.post.hmass.lg <- .brm.hmorph.lg %>%
  # extract posterior distribuitions
  posterior_samples() %>%
  # select only the main effects (start with "b_")
  select(grep("^b_hatchmass", colnames(.))) %>%
  # calculate the predicted posterior distributions for each treatment level
  reframe(Trait = "Hatchling Mass", # add trait column
          Species = "L. guichenoti", # add species column
          `Control,23°C` = b_hatchmass_Intercept, # Control, 23 is the default
          `Control,28°C` = b_hatchmass_Intercept + b_hatchmass_temp28, # add effect of temp
          `CORT,23°C` = b_hatchmass_Intercept + b_hatchmass_cortB, # add effect of cort
          `CORT,28°C` = b_hatchmass_Intercept + b_hatchmass_cortB + b_hatchmass_temp28) # add effect of temp and cort
#----
# -- hatchling svl 
#----
# -| best model: hatchling svl ~ cort + temp + egg mass + (1|clutch)
# --- lampropholis delicata
.post.hsvl.ld <- .brm.hmorph.ld %>%
  # extract posterior distributions
  posterior_samples() %>%
  # select only main effects (start with "b_")
  select(grep("^b_hatchsvl", colnames(.))) %>%
  # calculate the predicted posterior distributions for each treatment level
  reframe(Trait = "Hatchling SVL", # add a trait column
          Species = "L. delicata", # add a species column
          `Control,23°C` = b_hatchsvl_Intercept, # control,23 is the default
          `Control,28°C` = b_hatchsvl_Intercept + b_hatchsvl_temp28, # add effect of temp
          `CORT,23°C` = b_hatchsvl_Intercept + b_hatchsvl_cortB, # add effect of cort
          `CORT,28°C` = b_hatchsvl_Intercept + b_hatchsvl_cortB + b_hatchsvl_temp28) # add effect of cort and temp
# --- lampropholis guichenoti
.post.hsvl.lg <- .brm.hmorph.lg %>%
  posterior_samples() %>%
  select(grep("^b_hatchsvl", colnames(.))) %>%
  reframe(Trait = "Hatchling SVL",
          Species = "L. guichenoti",
          `Control,23°C` = b_hatchsvl_Intercept,
          `Control,28°C` = b_hatchsvl_Intercept + b_hatchsvl_temp28,
          `CORT,23°C` = b_hatchsvl_Intercept + b_hatchsvl_cortB,
          `CORT,28°C` = b_hatchsvl_Intercept + b_hatchsvl_cortB + b_hatchsvl_temp28)
#----
# -- merge extracted posteriors
.post <- bind_rows(.post.eggox.ld,
                   .post.eggox.lg,
                   .post.hmass.ld,
                   .post.hmass.lg)
# - run pairwise comparisons between each treatment group for each response variable
.tbl.s.2 <- .post %>%
  group_by(Trait, Species) %>%
  reframe(`[CORT,28°C] - [CORT,23°C]` = `CORT,28°C` - `CORT,23°C`,
          `[CORT,28°C] - [Control,28°C]` = `CORT,28°C` - `Control,28°C`,
          `[CORT,28°C] - [Control,23°C]` = `CORT,28°C` - `Control,23°C`,
          `[CORT,23°C] - [Control,28°C]` = `CORT,23°C` - `Control,28°C`,
          `[CORT,23°C] - [Control,23°C]` = `CORT,23°C` - `Control,23°C`,
          `[Control,28°C] - [Control,23°C]` = `Control,28°C` - `Control,23°C`
  ) %>%
  pivot_longer(cols = 3:ncol(.),
               names_to = "Contrast") %>%
  mutate(Contrast = factor(Contrast, levels = c("[CORT,28°C] - [CORT,23°C]",
                                                "[CORT,28°C] - [Control,28°C]",
                                                "[CORT,28°C] - [Control,23°C]",
                                                "[CORT,23°C] - [Control,28°C]",
                                                "[CORT,23°C] - [Control,23°C]",
                                                "[Control,28°C] - [Control,23°C]"))) %>%
  arrange(Contrast) %>%
  group_by(Trait, Species, Contrast) %>%
  # calculate model summary stats
  reframe(Estimate = mean(value), # mean of distribution
          Est.Error = se(value), # error about the mean of distribution
          `l-95% CI` = quantile(value, 0.025), # lower 95% CI
          `u-95% CI` = quantile(value, 0.975), # upper 95% CI
          pMCMC = pmcmc(value)) %>% # pMCMC value
  # round everything to 3 decimal places for easy viewing
  mutate(across(.cols = Estimate:pMCMC,
                .fns = ~round(.x, digits = 3)))

# split the supplementary table into components for summary RMD
.tbl.s.2.a <- .tbl.s.2 %>%
  filter(Trait == "Egg Mortality")
.tbl.s.2.a.i <- .tbl.s.2.a %>%
  filter(Species == "L. delicata") %>%
  .fx.ft()
.tbl.s.2.a.ii <- .tbl.s.2.a %>%
  filter(Species == "L. guichenoti") %>%
  .fx.ft()
.tbl.s.2.b <- .tbl.s.2 %>%
  filter(Trait == "Hatchling Mass")
.tbl.s.2.b.i <- .tbl.s.2.b %>%
  filter(Species == "L. delicata") %>%
  .fx.ft()
.tbl.s.2.b.ii <- .tbl.s.2.b %>%
  filter(Species == "L. guichenoti") %>%
  .fx.ft()
.tbl.s.2.c <- .tbl.s.2 %>%
  filter(Trait == "Hatchling SVL")
.tbl.s.2.c.i <- .tbl.s.2.c %>%
  filter(Species == "L. delicata") %>%
  .fx.ft()
.tbl.s.2.c.ii <- .tbl.s.2.c %>%
  filter(Species == "L. guichenoti") %>%
  .fx.ft()
# convert everything to flextables
.tbl.s.2.a <- .fx.ft(.tbl.s.2.a)
.tbl.s.2.b <- .fx.ft(.tbl.s.2.b)
.tbl.s.2.c <- .fx.ft(.tbl.s.2.c)
.tbl.s.2 <- .fx.ft(.tbl.s.2)
# save the supplementary table 2 (pairwise comparisons within species) flextable
save_as_image(.tbl.s.2, here("output", "tables", "tbl s2 - treatment pairwise comparisons.png"))
trt.pairwise <- .tbl.s.2
#----
# - run pairwise comparisons between species within each treatment group for each response variable
.tbl.3 <- .post %>%
  pivot_longer(3:ncol(.),
               names_to = "Treatment",
               values_to = "value") %>%
  group_by(Trait, Treatment, Species) %>%
  mutate(tracer = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = "Species",
              values_from = "value") %>%
  select(-tracer) %>%
  group_by(Trait, Treatment) %>%
  reframe(value = `L. delicata` - `L. guichenoti`) %>%
  group_by(Trait, Treatment) %>%
  reframe(Contrast = "[L. delicata] - [L. guichenoti]",
          Estimate = mean(value),
          Est.Error = se(value),
          `l-95% CI` = quantile(value, 0.025),
          `u-95% CI` = quantile(value, 0.975),
          pMCMC = pmcmc(value)) %>%
  mutate(across(.cols = Estimate:pMCMC,
                .fns = ~round(.x, digits = 3))) %>%
  .fx.ft()

# save the flextable for supplementary table 3 (species level pairwise comparisons)
save_as_image(.tbl.3, here("output", "tables", "tbl s3 - species comparisons.png"))

# checking differences in main effects between species
.tbl.3 <- .post %>%
  group_by(Trait, Species) %>%
  reframe(CORT = c(`CORT,23°C`, `CORT,28°C`),
          Control = c(`Control,23°C`, `Control,28°C`),
          `23°C` = c(`CORT,23°C`, `Control,23°C`),
          `28°C` = c(`CORT,28°C`, `Control,28°C`)) %>%
  group_by(Trait, Species) %>%
  reframe(`Stress\n[CORT - Control]` = CORT - Control,
          `Temperature\n[28°C - 23°C]` = `28°C` - `23°C`) %>%
  pivot_longer(cols = c(`Stress\n[CORT - Control]`, `Temperature\n[28°C - 23°C]`),
               names_to = "Effect",
               values_to = "value") %>%
  group_by(Trait, Species, Effect) %>%
  mutate(tracer = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = "Species",
              values_from = "value") %>%
  group_by(Trait, Effect) %>%
  reframe(Contrast = "L. delicata -\n L. guichenoti",
          value = `L. delicata` - `L. guichenoti`) %>%
  group_by(Trait, Effect, Contrast) %>%
  reframe(Estimate = mean(value),
          Est.Error = se(value),
          `l-95% CI` = quantile(value, 0.025),
          `u-95% CI` = quantile(value, 0.975),
          pMCMC = pmcmc(value)) %>%
  mutate(across(.cols = Estimate:pMCMC,
                .fns = ~round(.x, digits = 3))) %>%
  filter(Trait != "Hatchling SVL") %>%
  .fx.ft()
# save the flextable for species-level comparisons of main effects
save_as_image(.tbl.3, here("output", "tables", "tbl3 - between species main effects comparison.png"))
spp.pairwise <- .tbl.3



