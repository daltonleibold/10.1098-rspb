#---
#title: "func"
#author: "Dalton C. Leibold"
#date: 1/11/2023
#---

pacman::p_load(tidyverse)

#' @title pMCMC Function
#' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0.
#' @param null A numeric value decsribing what the null hypothesis should be
#' @param twotail Whether to conduct a one-tailed hypothesis or a two-tailed hypotheses. Default = true indicating a two-tailed test will be done.
pmcmc <- function(x, null = 0, twotail = TRUE){
  if(twotail){
    2*(1 - max(table(x<=null) / length(x)))
  } else{
    (1 - max(table(x<=null) / length(x)))
  }
}

#' @title re.code
#' @param x a character or factor vector
#' @param mapping a vector of matched channel names to label names
#' @description function for re-labelling channels with new parameter names
#' @export
re.code <- function(x, mapping) {
  labels <- base::unname(mapping)
  base::names(labels) <- base::names(mapping)
  out <- labels[base::match(x, base::names(mapping))]
  out[base::is.na(out)] <- x[base::is.na(out)]
  return(out)
}

# se = sd / sqrt(n)
# sd: population standard deviation
# n: number of elements in the population
se <- function(x, rm_na = TRUE){
  sd(x, na.rm = rm_na) / sqrt(length(x))
}

# making a function for quickly summarising model outputs
mod.sum <- function(mod){
  # get only the fixed effects
  summary(mod)[["fixed"]] %>%
    # make a pmcmc columns based on the posterior samples
    mutate(pMCMC = posterior_samples(mod) %>%
             # select only the fixed effects columns, which start with b_
             select(grep("^b_", colnames(.))) %>%
             # calculate the pmcmc value of each fixed effect
             # note that this only works if there are no interaction terms
             reframe(across(.cols = everything(),
                            .fns = pmcmc)) %>%
             # pivot to long form so we can extract the vector
             pivot_longer(cols = everything(),
                          names_to = "effect",
                          values_to = "pmcmc") %>%
             # extract the pmcmc vector
             # will be in same order as the model summary
             .$pmcmc) %>%
    # select only the columns we need to display
    select(Estimate, Est.Error, `l-95% CI`, `u-95% CI`, pMCMC) %>%
    # round to significant figures
    mutate(across(.cols = everything(),
                  .fns = ~round(.x, digits = 3))) %>%
    # suppress warnings because posterior_samples is deprecated
    suppressWarnings()
}


