-----
# - denotes a folder; multiple # (e.g., ####) indicate subdirectories within that folder
> - denotes a single file
$ - denotes a column within a csv or spreadsheet
----
# RSPB-2025-1226 - Repository for data, code, and supplementary materials for doi: 10.1098/rspb ('Life-history strategy mediates the effects of multiple developmental stressors on Australian skinks').

> `README.md` - a read me file explaining what is in this repo; you're reading it right now! 
> `refs.bib` - a bibtex file containing bibliographic information on all the papers cited; has been hand checked for rendering a references list but may still have some issues in the metadata

## data - folder containing raw data and spreadsheets of processed data that can be used to make tables
> data.csv - csv spreadsheet of raw data used in all downstream analyses
$species - a factor of which species the subject is; either 'delicata' or 'guichenoti'
$egg_id - a unique identification number for the individual animal
$clutch - an identifier of the clutch to which an egg belonged
$cort - a factor indicating the developmental stress treatment the individual received; either A (Control) or B (CORT); DCL, AYP, and PRS were blind to the treatments, hence the A and B notation
$temp - a factor indicating the incubation temperature treatment the individual received. either 23 (23 ± 3°C) or 28 (28 ± 3°C)
$egg_date - the date (dd/mm/yyyy) an egg was found
$egg_mass - the mass (in grams) of an egg upon being found
$egg_length - the length (in millimetres) along the longest portion of the egg upon being found
$egg_width - the length (in millimetres) along the widest portion of the egg upon being found
$hatch_date -  the date (dd/mm/yyyy) an egg hatched
$hatch_mass - the mass (in grams) of a hatchling upon emerging from an egg
$hatch_svl - the snout-vent-length (in millimetres) of a hatchling upon emerging from an egg
$hatch_tl - the tail-length (in millimetres) of a hatchling upon merging from an egg; not used in any analyses
$death_date - the date (dd/mm/yyyy) an individual died, if they died
$egg_death - a factor indicating whether an egg hatched ('hatched') or died prior to hatching ('died')
$egg_time - the total length of time (in days) an individual remained in the egg; if an egg died, this is the number of days until they died; if an egg hatched, this is the number of days until they hatched
$hatch_death - a factor indicating whether an individual remained alive ('alive') or died ('died') within the first few months following hatching; this was unused
$death_age - the age (in days) since hatching until an individual died if they died; this was unused
$egg_ox - a logical indicating whether an egg died (1) or hatched (0); mostly synonymous with egg_death, but was slightly more robust for analyses

## output - folder containing r script output like figures, tables, and models
### figures - folder containing figures; pretty self-explanatory
> `fig1 - methods schematic.png` - figure 1; see caption in manuscript for more details
> `fig2 - egg mortality and hatch mass.png` - figure 2; see caption in manuscript for more details
#### `fig2 - data plots` - folder containing components to make figure 2 and an expanded version of figure 2
> `fig_2_a_i - model egg mortality.png` - fadecloud plot of egg mortality values by treatment and species; data is posterior distributions from model `brm_eggox_ld.rds` and `brm_eggox_lg.rds`
> `fig_2_a_ii - raw egg mortality.png` - bar chart of egg mortality by treatment and species; data is raw counts of eggs that either died or hatched
> `fig_2_b_i - model hatchling mass.png` - fadecloud plot of hatchling mass values by treatment and species; data is predicted posterior distributions from model `brm_hmorph_ld.rds` and `brm_hmorph_lg.rds`
> `fig_2_b_ii - raw hatchling mass.png` - summary statistic (mean and se) plot of hatchling mass values by treatment and species; data is from raw measurements of hatchling mass
> `fig_2_c_i - model hatchling svl.png` - fadecloud plot of hatchling svl values by treatment and species; data is predicted posterior distributions from model `brm_hmorph_ld.rds` and `brm_hmorph_lg.rds`
> `fig_2_c_ii - raw hatchling svl.png` - summary statistic (mean and se) plot of hatchling svl values by treatment and species; data is from raw measurements of hatchling svl
### models - folder containing models
> `eggox_ld.rds` - bayesian regression generalized linear mixed-effects model model of l. delicata egg mortality as a function of temp, cort, and egg mass with clutch as a random factor
> `eggox_lg.rds` - bayesian regression generalized linear mixed-effects model model of l. guichenoti egg mortality as a function of temp, cort, and egg mass with clutch as a random factor
> `hmorph_ld.rds` - multivariate bayesian regression generalized linear mixed-effects model of l. delicata hatchling mass and hatchling svl as a function of temp, cort, and egg mass with clutch as a random factor
 > `hmorph_lg.rds` - multivariate bayesian regression generalized linear mixed-effects model of l. guichenoti hatchling mass and hatchling svl as a function of temp, cort, and egg mass with clutch as a random factor
### tables - folder containing tables; pretty self-explanatory
#### `tbl_s4 - model comparisons` - folder containing model comparison flextables used in creating `supplementary document - model comparisons.docx`
> `tbl s1 - all summary stats.png` - supplementary table 1; see `Supplementary Document - Supplementary Figures and Captions.docx`
> `tbl s2 - treatment pairwise comparisons.png` - supplementary table 2; see `Supplementary Document - Supplementary Figures and Captions.docx`
> `tbl s3 - species pariwise comparisons.png` - supplementary table 3; see `Supplementary Document - Supplementary Figures and Captions.docx`
> `tbl1 - summary stats.png` - table 1; see manuscript
> `tbl2 - stats test.png` - table 2; see manuscript
> `tbl3 - between species main effects comparison.png` - table 3; see manuscript

## r - folder containing r scripts for processing, analyzing, visualizing, and summarizing data
> func.R - r script containing custom functions used in our analyses. includes code snippets for calculating pmcmc values (our metric of significance in stats tests), calculating standard error, relabelling factor columns,
> hypothesis_tests.R - r script for performing hypothesis tests of main effects on embryonic mortality, hatchling mass, and hatchling svl in lampropholis delicata and lampropholis guichenoti. creates table 2 and table 3.
> model_comparisons.R - r script for performing model comparisons to decide the models of best fit for analyses. it then runs the models of best fit, which are stored in the 'models' folder. creates all the figures needed for 'supplementary document - model comparisons'. 
> plot_making.R - r script for making raw figures from data. creates fadecloud plots of egg mortality, hatchling mass, and hatchling svl divided by each treatment and each species. also creates summary stats plots for hatchling mass and hatchling svl, and a bar chart for proportions of egg mortality in each treatment group for each species. creates components of figure 2.
> summary_stats.R - r script for calculating summary statistics from the raw data. creates table 1 and supplementary table 1.
> source.R - r script that processes the data and sources in all of the relevant information from hypothesis_tests.R, plot_making.R, and summary_stats.R




