# summary.R
# script for summarising raw data
# prepares tbl1

#  load libraries
pacman::p_load(tidyverse, flextable)
# set working directory to project root
setwd(here::here())
# read data
data <- read_csv(here::here("data", "data.csv"))
# source custom functions
source(here::here("r", "func.R"))
# set flextable defaults
set_flextable_defaults(
  font.family = "Arial",
  font.size = 12,
  digits = 3
)
# make a basic function for making flextables
.fx.ft <- function(x){
  x %>%
    flextable() %>%
    padding(padding = 10, part = "all") %>%
    width(width = 1.5) %>%
    width(width = 2.5,
          j = 2) %>%
    merge_v(j = 1:2) %>%
    valign(j = 1:2, valign = "top") %>%
    align(align = "center", part = "all") %>%
    bg(bg = "white", part = "all")
}
.tbl.1.a <- data %>%
  # remove eggs we are unsure of
  filter(!is.na(egg_death)) %>%
  # make counts of the numbers of eggs
  group_by(species, cort, temp, egg_death) %>%
  reframe(n = n()) %>%
  # pivot to wide form
  pivot_wider(names_from = "egg_death",
              values_from = "n") %>%
  # calculate the proportion of eggs that died
  mutate(total = died + hatched,
         egg_mortality = died/total * 100) %>%
  # can drop the total column now
  select(-total) %>%
  # give them more informative names
  rename(`eggs (died)` = died,
         `eggs (hatched)` = hatched) %>%
  # pivot back to long form
  pivot_longer(`eggs (died)`:egg_mortality,
               names_to = "Trait",
               values_to = "value") %>%
  # format trait names
  mutate(Trait = str_replace(Trait, "_", " "),
         Trait = str_to_title(Trait)) %>%
  # round the percentages
  mutate(value = ifelse(Trait == "Egg Mortality",
                        round(value, digits = 1),
                        round(value, digits = 0)),
         # add the units (%)
         value = ifelse(Trait == "Egg Mortality",
                        str_c(value, "%", sep = ""),
                        value)) %>%
  # format treatment and species names
  mutate(cort = re.code(cort, c("A" = "Control", "B" = "CORT")),
         temp = re.code(temp, c("23" = "23°C", "28" = "28°C")),
         species = str_c("L. ", species, sep = "")) %>%
  # join cort and temperature treatments
  unite("Treatment", c(cort, temp), sep = ",") %>%
  # pivot to wide form by treatment combination
  pivot_wider(names_from = "Treatment",
              values_from = "value") %>%
  # capitalise species
  rename(Species = species) %>%
  # change everything to character for joining tables
  mutate(across(.cols = everything(),
                .fns = as.character))

.tbl.1.b <- data %>%
  # remove eggs we are unsure of
  filter(!is.na(egg_death)) %>%
  # calculate mean and standard deviation of egg morphology traits
  group_by(species, cort, temp) %>%
  reframe(across(.cols = c(egg_mass:egg_width),
                 .fns = list(mean = ~mean(.x, na.rm = T),
                             sd = ~sd(.x, na.rm = T)))) %>%
  # pivot and separate to tabulate
  pivot_longer(egg_mass_mean:egg_width_sd) %>%
  separate(name, into = c("stage", "trait", "stat"), sep = "_") %>%
  unite(trait, c("stage", "trait"), sep = " ") %>%
  # re-code columns
  mutate(trait = str_to_title(trait),
         species = str_c("L.", species, sep = " "),
         cort = re.code(cort, c("A" = "Control", "B" = "CORT")),
         temp = re.code(temp, c("23" = "23°C", "28" = "28°C"))) %>%
  unite(treatment, c(cort, temp), sep = ",") %>%
  # round values
  mutate(value = round(value, digits = 3)) %>%
  # pivot back to wide form
  pivot_wider(names_from = "stat",
              values_from = "value") %>%
  # unite mean and standard deviations
  unite("value", c(mean, sd), sep = "\n± ") %>%
  # pivot even wider by treatment level
  pivot_wider(names_from = "treatment",
              values_from = "value") %>%
  # format column names
  rename(Species = species,
         Trait = trait) %>%
  # change to character for joining
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  # specify that this includes all eggs
  mutate(Trait = str_c(Trait, "(All)", sep = " "))

.tbl.1.c <- data %>%
  # filter out eggs we are unsure of
  filter(!is.na(egg_death)) %>%
  # group by whether eggs died or not this time
  group_by(species, cort, temp, egg_death) %>%
  # calculate the mean and standard deviation of each egg morphology trait
  reframe(across(.cols = c(egg_mass:egg_width),
                 .fns = list(mean = ~mean(.x, na.rm = T),
                             sd = ~sd(.x, na.rm = T)))) %>%
  # pivot to long form
  pivot_longer(egg_mass_mean:egg_width_sd) %>%
  # split into stage, trait, and statistic
  separate(name, into = c("stage", "trait", "stat"), sep = "_") %>%
  unite(trait, c("stage", "trait"), sep = " ") %>%
  # re-code columns to publication naming convention
  mutate(trait = str_to_title(trait),
         egg_death = str_to_title(egg_death),
         egg_death = str_c("(", egg_death, ")", sep = ""),
         species = str_c("L.", species, sep = " "),
         cort = re.code(cort, c("A" = "Control", "B" = "CORT")),
         temp = re.code(temp, c("23" = "23°C", "28" = "28°C"))) %>%
  # re-format names and collapse columns
  unite(trait, c(trait, egg_death), sep = " ") %>%
  unite(treatment, c(cort, temp), sep = ",") %>%
  # round values for easier viewing
  mutate(value = round(value, digits = 3)) %>%
  # pivot back to wide form
  pivot_wider(names_from = "stat",
              values_from = "value") %>%
  # join the mean and standard deviation
  unite("value", c(mean, sd), sep = "\n± ") %>%
  # pivot even wider by the treatment combination
  pivot_wider(names_from = "treatment",
              values_from = "value") %>%
  # reformat column names
  rename(Species = species,
         Trait = trait) %>%
  # convert everything to character for joining
  mutate(across(.cols = everything(),
                .fns = as.character))

.tbl.1.d <- data %>%
  # calculate mean and standard deviation of each hatchling morphology trait
  group_by(species, cort, temp) %>%
  reframe(across(.cols = hatch_mass:hatch_tl,
                 .fns = list(mean = ~mean(.x, na.rm = T),
                             sd = ~sd(.x, na.rm = T)))) %>%
  # pivot to long form
  pivot_longer(hatch_mass_mean:hatch_tl_sd) %>%
  # split columns into stage, trait, and statistic
  separate(name, into = c("stage", "trait", "stat"), sep = "_") %>%
  # round values for easier viewing
  mutate(value = round(value, digits = 3)) %>%
  # pivot to wide form
  pivot_wider(names_from = "stat",
              values_from = "value") %>%
  # join the mean and standard deviation
  unite("value", c(mean, sd), sep = "\n± ") %>%
  # recode columns
  mutate(cort = re.code(cort, c("A" = "Control", "B" = "CORT")),
         temp = re.code(temp, c("23" = "23°C", "28" = "28°C"))) %>%
  mutate(stage = str_to_title(stage),
         trait = ifelse(trait == "mass", "Mass", str_to_upper(trait))) %>%
  # join the treatment columns
  unite("Treatment", c(cort, temp), sep = ",") %>%
  # unite the stage and trait columns
  unite("Trait", c(stage, trait), sep = " ") %>%
  # pivot to wide form by treatment combination
  pivot_wider(names_from = "Treatment",
              values_from = "value") %>%
  # reformat species name
  mutate(species = str_c("L. ", species, sep = "")) %>%
  rename(Species = species) %>%
  # add units
  mutate(Trait = ifelse(Trait == "Hatch Mass",
                        str_c(Trait, "(g)", sep = " "),
                        str_c(Trait, "(mm)", sep = " "))) %>%
  # convert everything to character for joining
  mutate(across(.cols = everything(),
                .fns = as.character))

# make a full table of all the summary statistics
.tbl.s.1 <- .tbl.1.a %>%
  # bind all the tables together
  bind_rows(.tbl.1.b, .tbl.1.c, .tbl.1.d) %>%
  # arrange by species
  arrange(Species) %>%
  # better names
  mutate(Trait = str_replace(Trait, "Hatch Mass", "Hatchling Mass"),
         Trait = str_replace(Trait, "Hatch SVL", "Hatchling SVL"),
         Trait = str_replace(Trait, "Hatch TL", "Hatchling TL")) %>%
  # rearrange columns
  select(Species, Trait, 
         `Control,23°C`, `CORT,23°C`, 
         `Control,28°C`, `CORT,28°C`)

# make a subset for tbl.1
.tbl.1 <- .tbl.s.1 %>%
  filter(Trait == "Eggs (Died)"
         | Trait == "Eggs (Hatched)"
         | Trait == "Egg Mortality"
         | Trait == "Egg Mass (Died)"
         | Trait == "Egg Mass (Hatched)"
         | Trait == "Hatchling Mass (g)"
         )

# convert each table to a flextable
.tbl.s.1 <- .fx.ft(.tbl.s.1)
.tbl.1 <- .fx.ft(.tbl.1)

# save the flextables
save_as_image(.tbl.s.1, here("output", "tables", "tbl s1 - all summary stats.png"))
save_as_image(.tbl.1, here("output", "tables", "tbl1 - summary stats.png"))

summary.stats <- .tbl.s.1
