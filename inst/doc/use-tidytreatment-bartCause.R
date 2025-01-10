## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(6, 4)
)

suppressPackageStartupMessages({
    library(bartCause)
    library(stan4bart)
    library(tidytreatment)
    library(dplyr)
    library(tidybayes)
    library(ggplot2)
  })
  
  # load pre-computed data and model
  sim <- suhillsim2_ranef

  

## ----load-data-print, echo = TRUE, eval = FALSE-------------------------------
# 
# # load packages
# library(bartCause)
# library(stan4bart)
# library(tidytreatment)
# library(dplyr)
# library(tidybayes)
# library(ggplot2)
# 
# # set seed so vignette is reproducible
# set.seed(101)
# 
# # simulate data
# sim <- simulate_su_hill_data(n = 100, treatment_linear = FALSE,  omega = 0, add_categorical = TRUE,
#                              n_subjects = 10, sd_subjects = 0.1,
#                              coef_categorical_treatment = c(0,0,1),
#                              coef_categorical_nontreatment = c(-1,0,-1)
# )
# 

## ----data-summary, echo = TRUE, eval = TRUE-----------------------------------

# non-treated vs treated counts:
table(sim$data$z)

dat <- sim$data
# a selection of data
dat %>% select(y, z, c1, x1:x3) %>% head()

# repeated observation counts for subjects:
table(sim$data$subject_id)


## ----run-bart, echo = TRUE, eval = TRUE---------------------------------------
  
# STEP 1 VS Model: Regress y ~ covariates
vs_bart <- stan4bart(y ~ bart(. - subject_id - z) + (1|subject_id), 
                             data = dat, iter = 5000, verbose = -1)

# STEP 2: Variable selection
  # select most important vars from y ~ covariates model
  # very simple selection mechanism. Should use cross-validation in practice
covar_ranking <- covariate_importance(vs_bart)
var_select <- covar_ranking %>% 
  filter(avg_inclusion > mean(avg_inclusion) - sd(avg_inclusion)) %>% # at minimum: within 1 sd of mean inclusion
  pull(variable)

# change categorical variables to just one variable
var_select <- unique(gsub("c1.[1-3]$","c1", var_select))

var_select
# includes all covariates

# STEP 3 PS Model: Regress z ~ selected covariates
ps_bart <- stan4bart(z ~ bart(. - subject_id - y) + (1|subject_id), 
                             data = dat, iter = 5000, verbose = -1)

# store propensity score in data
prop_score <- fitted(ps_bart)

# Step 4 TE Model: Regress y ~ z + covariates + propensity score
te_bart <- bartc(response = y, treatment = z, confounders = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,  
                 parametric = (1|subject_id), data = dat, method.trt = prop_score, 
                 iter = 5000, bart_args = list(keepTrees = TRUE))

#* The posterior samples are kept small to manage size on CRAN


## ----tidy-bart-fit, echo=TRUE, cache=FALSE------------------------------------

# get model parameters (excluding BART paramaters)
posterior_params <- tidy_draws(te_bart)

posterior_fitted <- epred_draws(te_bart, value = "fitted")


## ----tidy-bart-pred, eval=FALSE, echo=TRUE, cache=FALSE-----------------------
# 
# # Function to tidy predicted draws (adds predicted noise to fitted values)
# posterior_pred <- predicted_draws(te_bart, value = "predicted")
# 

## ----plot-tidy-bart, echo=TRUE, cache=FALSE-----------------------------------

treatment_var_and_c1 <- 
  dat %>% 
  select(z,c1) %>%
  mutate(.row = 1:n(), z = as.factor(z))

posterior_fitted %>%
  left_join(treatment_var_and_c1, by = ".row") %>%
  ggplot() + 
  stat_halfeye(aes(x = z, y = fitted)) + 
  facet_wrap(~c1, labeller = as_labeller( function(x) paste("c1 =",x) ) ) +
  xlab("Treatment (z)") + ylab("Posterior predicted value") +
  theme_bw() + ggtitle("Effect of treatment with 'c1' on posterior fitted values")


## ----post-treatment, eval = T-------------------------------------------------

# sample based (using data from fit) conditional treatment effects, posterior draws
posterior_treat_eff <- treatment_effects(te_bart)

# check lines up with summary results...


## ----cates-hist, echo=TRUE, cache=FALSE, eval = T-----------------------------

# Histogram of treatment effect (all draws)
posterior_treat_eff %>% 
  ggplot() +
  geom_histogram(aes(x = icate), binwidth = 0.1, colour = "white") + 
  theme_bw() + ggtitle("Histogram of treatment effect (all draws)")

# Histogram of treatment effect (median for each subject)
posterior_treat_eff %>% summarise(cte_hat = median(icate)) %>%
  ggplot() +
  geom_histogram(aes(x = cte_hat), binwidth = 0.1, colour = "white") + 
  theme_bw() + ggtitle("Histogram of treatment effect (median for each subject)")


