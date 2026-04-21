conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(brms::ar)
conflicted::conflicts_prefer(dplyr::lag)
conflicted::conflicts_prefer(kableExtra::group_rows)

##########################################################################################
##                     Study and selection of the different growth curves               ##
##                                   (Moulis version)                                   ##
##                             "Genetics of growth analysis"                            ##
##                          Anaïs Aragon & Pierre de Villemereuil                       ##
##                                        (2023)                                        ##
##########################################################################################

## Packages required
library(tidyverse)
library(here)
library(brms)
library(future)
library(readODS)
library(kableExtra)

## Source the custom functions for this analysis
source("_functions.R")

## Setting up a (lucky) seed
seed <- 777
set.seed(seed)

## Setting up the ggplot2 theme
theme_set(theme_bw())
# Colours for males and females
col_sex <- c("#55aa00", "#00007f")
col_sex_points <- c("#a4ff7d44", "#5560ff44")

## Setting up parallelisation
ncores <- min(parallel::detectCores() - 2, 20)
options(mc.cores = ncores)
plan(multicore)

## Setting up the number of chains
# For unique model to compute LOOIC
nchains <- 10

## Choosing parts of the script to execute
run_models_selection <- FALSE


## --------------------------------------- Loading and formatting the data ----

# Loading the Capture dataset
tbl_capt <- readRDS(here("Data/Moulis_capture_clean.rds"))

# Filtering the data:
# - Removing missing Sex values
tbl_capt <-
    tbl_capt |>
    filter(!is.na(Sex)) |>
    filter(n_distinct(Age_Year) >= 2, .by = ID)

## ------------------------------------------ Running models for selection ----

## ---- Logistic growth model ----

# Priors for the fixed effects of the non-linear model
# (Note that L0 is set as the size at birth)
# We used prior knowledge in the bibliography to setup informative priors
prior_log <-
    prior(normal(60, 10), nlpar = "Lmax") +
    prior(normal(0, 10), nlpar = "Lmax", coef = "SexM") +
    prior(uniform(0, 10), nlpar = "k") +
    prior(normal(0, 5), nlpar = "k", coef = "SexM")

# Setting up the formula
form_log <-
    bf(Size ~ Lmax  / (1 + ((Lmax - Size_NB) / Size_NB) * exp(- k * Age)),
    Lmax + k ~ Sex + (1|ID),
    nl = TRUE)

# Running the brms model
if (run_models_selection) {
    print("Running the logistic growth model")
    mod_log <-
        brm(formula = form_log,
            data    = tbl_capt,
            chains  = nchains,
            cores   = ncores,
            init    = rep(list(list(b_Lmax    = array(data = c(68, 0)),
                                    b_k       = array(data = c(1, 0)))),
                            nchains),
            prior   = prior_log,
            seed    = seed,
            iter    = 3000,
            warmup  = 1000,
            thin    = 5)
    saveRDS(mod_log, here("Output/Moulis_comparison_LOO_model_log.rds"), compress = "xz")
}

## ---- von Bertalanffy model ----

# Priors for the fixed effects of the non-linear model
# (Note that L0 is set as the size at birth)
# We used prior knowledge in the bibliography to setup informative priors
prior_vb <-
    prior(normal(60, 10), nlpar = "Lmax") +
    prior(normal(0, 10), nlpar = "Lmax", coef = "SexM") +
    prior(uniform(0, 10), nlpar = "k") +
    prior(normal(0, 5), nlpar = "k", coef = "SexM")

# Setting up the formula
form_vb <-
    bf(Size ~ Lmax - (Lmax - Size_NB) * exp(- k * Age),
    Lmax + k ~ Sex + (1|ID),
    nl = TRUE)

# Running the brms model
if (run_models_selection) {
    print("Running the von Bertalanffy growth model")
    mod_vb <-
        brm(formula = form_vb,
            data    = tbl_capt,
            chains  = nchains,
            cores   = ncores,
            init    = rep(list(list(b_Lmax    = array(data = c(68, 0)),
                                    b_k       = array(data = c(1, 0)))),
                          nchains),
            prior   = prior_vb,
            seed    = seed,
            iter    = 3000,
            warmup  = 1000,
            thin    = 5)
    saveRDS(mod_vb, here("Output/Moulis_comparison_LOO_model_vb.rds"), compress = "xz")
}

## ---- Gompertz model ----

# Priors for the fixed effects of the non-linear model
# (Note that L0 is set as the size at birth)
# We used prior knowledge in the bibliography to setup informative priors
prior_gz <-
    prior(normal(60, 10), nlpar = "Lmax") +
    prior(normal(0, 10), nlpar = "Lmax", coef = "SexM") +
    prior(uniform(0, 10), nlpar = "k") +
    prior(normal(0, 5), nlpar = "k", coef = "SexM")

# Setting up the formula
form_gz <-
    bf(Size ~ Lmax * exp(log(Size_NB / Lmax) * exp(- k * Age)),
    Lmax + k ~ Sex + (1|ID),
    nl = TRUE)

# Running the brms model
if (run_models_selection) {
    print("Running the Gompertz growth model")
    mod_gz <-
        brm(formula = form_gz,
            data    = tbl_capt,
            chains  = nchains,
            cores   = ncores,
            init    = rep(list(list(b_Lmax    = array(data = c(68, 0)),
                                    b_k       = array(data = c(1, 0)))),
                          nchains),
            prior   = prior_gz,
            seed    = seed,
            iter    = 3000,
            warmup  = 1000,
            thin    = 5)
    saveRDS(mod_gz, here("Output/Moulis_comparison_LOO_model_gz.rds"), compress = "xz")
}


## ---------------------------- Model comparison using Leave-One-Out (LOO) ----

mod_log <- readRDS(here("Output/Moulis_comparison_LOO_model_log.rds"))
mod_vb  <- readRDS(here("Output/Moulis_comparison_LOO_model_vb.rds"))
mod_gz  <- readRDS(here("Output/Moulis_comparison_LOO_model_gz.rds"))

# # Using the LOOIC from the loo() function
# loo_comp <- loo(mod_log, mod_vb, mod_gz)

# Using the WAIC to avoid issues with PSIS
waic_comp   <- waic(mod_log, mod_vb, mod_gz)

if (rownames(waic_comp[["diffs"]])[1] != "mod_log") {
    stop("The logistic growth is not the best model anymore,\r
         you need to manually account for model selection.")
}

# saveRDS(loo_comp, here("Output/Moulis_comparison_LOO.rds"), compress = TRUE)
saveRDS(waic_comp, here("Output/Moulis_comparison_WAIC.rds"), compress = TRUE)

# Outputting a nice table
out_waic <-
    waic_comp[["diffs"]] |>
    as.data.frame() |>
    rownames_to_column(var = "Model") |>
    select(Model,
           WAIC         = waic,
           Delta_WAIC   = elpd_diff,
           SE           = se_diff) |>
    mutate(Model = recode(Model,
                          mod_log = "Logistic",
                          mod_gz  = "Gompertz",
                          mod_vb  = "Von Bertanlanffy"))

write_ods(out_waic,
          path = here("Tables/Moulis_WAIC.ods"))

## ------------------------------------ Graphical comparison of all models ----

# Setting up a dodge width for everything
dodge_width <- 0.5

# Getting the average size at birth for females and males
mean_nb_size <-
    tbl_capt |>
    filter(Age == 0) |>
    summarise(Size = mean(Size), .by = Sex) |>
    deframe()

# Generating the "background" of the data for the plot
p_base <-
    ggplot() +
    geom_count(aes(x = Age,
                   y = Size,
                   fill = Sex,
                   colour = Sex),
               shape = 21,
#                alpha = 0.2,
#                fill = "#FFFFFF00",
               position = position_dodge(width = dodge_width),
               data = tbl_capt) +
    scale_x_continuous(breaks = scales::breaks_pretty(),
                       limits = c(-0.5,10.2)) +
    scale_size_area(max_size = 3,
                    transform = "sqrt",
                    limits = c(1, 200),
                    breaks = c(1, 5, 20, 50, 100, 200)) +
    scale_fill_manual(values = col_sex_points) +
    scale_colour_manual(values = col_sex) +
    ylim(c(15, 80))

# A function to generate the parametrised growth curves
generate_growth_functions <- function(func_, sex_, l0_, lmax_, k_) {

    # Create a call to the func_ function with the parameter values
    body <- parse(text = str_c(func_,
                                    "(x, L0 = ", l0_,
                                    ", Lmax = ", lmax_,
                                    ", k = ", k_,
                                    ")"))

    # Generating the list of arguments without defaults (bit of a hack)
    args <- list()
    args[["x"]] <- alist(x=)$x

    # Generate the function and return it
    fun <- eval(call("function", as.pairlist(args), body[[1]]), parent.frame())
    return(fun)
}

# Getting the parameters from the models
tbl_func_comp <-
    crossing(
        tibble(Model_Name   = c("Logistic", "Bertalanffy", "Gompertz"),
               Function     = c("f_log", "f_vb", "f_gz")),
        tibble(Sex          = c("Female", "Male"),
               L0           = mean_nb_size[c("F", "M")],
               Colour       = col_sex)
    ) |>
    mutate(Model_Fit = case_when(Model_Name == "Logistic"    ~ list(mod_log),
                                 Model_Name == "Bertalanffy" ~ list(mod_vb),
                                 Model_Name == "Gompertz"    ~ list(mod_gz)),
           Fixed = map2_dfr(Model_Fit, Sex,
                            \(fit, sex) {
                                fix <- fixef(fit)
                                base <- c(Lmax  = fix["Lmax_Intercept", "Estimate"],
                                          k     = fix["k_Intercept", "Estimate"])
                                out <- base
                                if (sex == "Male") {
                                    out <- base +
                                           c(Lmax = fix["Lmax_SexM", "Estimate"],
                                             k    = fix["k_SexM", "Estimate"])
                                }
                                out <- as.list(out) |> as_tibble()
                                return(out)
                            })) |>
    unpack(Fixed) |>
    select(-Model_Fit) |>
    mutate(Function = pmap(list(func_   = Function,
                                sex_    = Sex,
                                l0_     = L0,
                                lmax_   = Lmax,
                                k_      = k),
                           generate_growth_functions),
           Plot = pmap(list(func_ = Function,
                            col_  = Colour,
                            name_ = Model_Name,
                            sex_  = Sex),
                       \(func_, col_, name_, sex_) {
                          geom_function(fun     = func_,
                                        colour  = col_,
                                        data    = tibble(Name = name_,
                                                         Sex  = sex_))
                       }))

p_comp <-
    p_base +
    tbl_func_comp[["Plot"]] +
    facet_grid(Name ~ .) +
    theme(legend.position = "bottom") +
    labs(y = "Size (mm)", x = "Age (years)")

cairo_pdf(here("Figures/Moulis_Comparison_Growth_Curves.pdf"), width = 7, height = 8)
plot(p_comp)
dev.off()

## ----------------------------- Outputting the graph for the final model ----

# Loading the fit of the final model
mod_final <- readRDS(here("Output/Moulis_comparison_LOO_model_log.rds"))

# Generating a small table for the male/female growth functions
tbl_func_final <-
    crossing(
        tibble(Model_Name   = "Logistic",
               Function     = "f_log"),
        tibble(Sex          = c("Female", "Male"),
               L0           = mean_nb_size[c("F", "M")],
               Colour       = col_sex)
    ) |>
    mutate(Model_Fit = list(mod_final, mod_final),
           Fixed = map2_dfr(Model_Fit, Sex,
                            \(fit, sex) {
                                fix <- fixef(fit)
                                base <- c(Lmax  = fix["Lmax_Intercept", "Estimate"],
                                          k     = fix["k_Intercept", "Estimate"])
                                out <- base
                                if (sex == "Male") {
                                    out <- base +
                                           c(Lmax = fix["Lmax_SexM", "Estimate"],
                                             k    = fix["k_SexM", "Estimate"])
                                }
                                out <- as.list(out) |> as_tibble()
                                return(out)
                            })) |>
    unpack(Fixed) |>
    select(-Model_Fit) |>
    mutate(Function = pmap(list(func_   = Function,
                                sex_    = Sex,
                                l0_     = L0,
                                lmax_   = Lmax,
                                k_      = k),
                           generate_growth_functions),
           Plot = pmap(list(func_ = Function,
                            col_  = Colour,
                            name_ = Model_Name,
                            sex_  = Sex),
                       \(func_, col_, name_, sex_) {
                          geom_function(fun     = func_,
                                        colour  = col_,
                                        data    = tibble(Name = name_,
                                                         Sex  = sex_))
                       }))

p_final <-
    p_base +
    tbl_func_final[["Plot"]] +
    theme(legend.position = "bottom") +
    labs(y = "Size (mm)", x = "Age (years)")

saveRDS(p_final, file = here("Output/Object_Final_Growth_Curve_Moulis.rds"))

cairo_pdf(here("Figures/Moulis_Final_Growth_Curve.pdf"), width = 7, height = 5)
plot(p_final)
dev.off()
