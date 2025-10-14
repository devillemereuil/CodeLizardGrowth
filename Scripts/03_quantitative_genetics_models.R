conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::lag)
conflicted::conflicts_prefer(dplyr::where)
conflicted::conflicts_prefer(brms::ar)
conflicted::conflicts_prefer(brms::me)
conflicted::conflicts_prefer(tidyr::expand)
conflicted::conflicts_prefer(tidyr::pack)
conflicted::conflicts_prefer(tidyr::unpack)
conflicted::conflicts_prefer(kableExtra::group_rows)

##########################################################################################
##                        Quantitative genetics of the growth curves                    ##
##                             "Genetics of growth analysis"                            ##
##                          Anaïs Aragon & Pierre de Villemereuil                       ##
##                                        (2023)                                        ##
##########################################################################################

## Packages required
library(tidyverse)
library(here)
library(brms)
library(MCMCglmm)
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

## Setting up parallelisation
n_cores  <- min(parallel::detectCores() - 2, 20)
options(mc.cores = n_cores)
plan(multicore)

## Choosing parts of the script to execute
run_models_brms <- FALSE
run_models_MCMCglmm <- FALSE

## Setting up the MCMC chains parameters for the linear models
n_chains_lm <- 10
n_iter_lm   <- 3000
n_wup_lm    <- 1000
thin_lm     <- 1

## Setting up the MCMC chains parameters for multiple imputation models
n_chains_nl <- 1
n_iter_nl   <- 3000
n_wup_nl    <- 1000
thin_nl     <- 1

## Increasing the size of globals to 2GiB for the multiple imputation model
options(future.globals.maxSize= 2 * 1024^3)

## --------------------------------------- Loading and formatting the data ----

## Capture dataset ----

# Loading the dataset
tbl_capt <- readRDS(here("Data/capture_avg_imputed.rds"))
lst_mult_capt <- readRDS(here("Data/capture_list_imputed.rds"))

# Generating a function to format the datasets
format_capt <- function(tbl_) {
    tbl_ <-
        tbl_ |>
        mutate(Age      = if_else(Stage == "JC", 0.2, Age),
               Year_fac = as.factor(Year)) |>
        filter(Stage != "JC")

    return(tbl_)
}
# Applying the formatting to all Capture datasets
tbl_capt <- format_capt(tbl_capt)
lst_mult_capt <- map(lst_mult_capt, format_capt)

## Pedigree dataset & relatedness matrix ----

tbl_ped <- readRDS(here("Data/pedigree_clean.rds"))
mat_A   <- readRDS(here("Data/relatedness_matrix.rds"))

## Reproduction dataset ----

tbl_repro <- readRDS(here("Data/repro_clean.rds"))
# Get size, mass and age of mothers at year of reproduction
tbl_repro <-
    tbl_repro |>
    left_join(tbl_capt |>
              select(ID_Mom = ID,
                     Year,
                     Age_Mom = Age,
                     Size_Mom = Size,
                     Weight_Mom = Weight) |>
              summarise(across(everything(), \(vec) mean(vec, na.rm = TRUE)),
                        .by = c(ID_Mom, Year)))

## ------------------------------------------ Model for Size at Birth (L0) ----

# Formatting a data for this analysis, esp. removing imputed size at birth
tbl_L0 <-
    tbl_capt |>
    filter(ID %in% rownames(mat_A),
           Stage == "NB",
           Size_SE == 0) |>
    left_join(tbl_repro |>
              select(ID = ID_Juv,
                     ID_Mom,
                     Date_Birth)) |>
    mutate(ID_Mom = if_else(is.na(ID_Mom), str_c("X_", 1:n()), ID_Mom),
           Clutch_Size = map2_dbl(ID_Mom, Year,
                                  \(id, yr) { sum(tbl_repro[["ID_Mom"]] == id &
                                                  tbl_repro[["Year"]] == yr) }),
           Date_Birth = str_c(year(today()),
                              month(Date_Birth),
                              day(Date_Birth),
                              sep = "-") |>
                        ymd(),
           Days_Birth = Date_Birth -
                        ymd(str_c(year(today()), "06", "01", sep = "-"))) |>
    mutate(Days_Birth = Days_Birth - mean(Days_Birth),
           .by = Year)

# Subsetting the relatedness matrix
mat_A_L0 <- mat_A[tbl_L0[["ID"]], tbl_L0[["ID"]]]

# Running the animal model for size at birth (L0)
form_L0 <- bf(Size ~ Sex + Days_Birth +
                     (1 | gr(ID, cov = A)) +
                     (1 | Year_fac) +
                     (1 | ID_Mom))

if (run_models_brms) {
    print("Running brms model for L0")

    mod_L0 <-
        brm(formula     = form_L0,
            data        = tbl_L0,
            data2       = list(A = as.matrix(mat_A_L0)),
            chains      = n_chains_lm,
            cores       = min(n_chains_lm, n_cores),
            seed        = seed,
            iter        = n_iter_lm,
            warmup      = n_wup_lm,
            thin        = thin_lm,
            save_pars   = save_pars(group = FALSE),
            control     = list(adapt_delta = 0.8))
    saveRDS(mod_L0, file = here("Output/animal_model_L0.rds"), compress = "xz")
}

if (run_models_MCMCglmm) {
    print("Running MCMCglmm model for L0")
    tbl_L0[["animal"]] <- tbl_L0[["ID"]]
    pior_Mg <- list(R = list(V = 1, nu = 0.02),
                    G = list(G1 = list(V = 1, nu = 0.02),
                             G2 = list(V = 1, nu = 0.02),
                             G3 = list(V = 1, nu = 0.02)))
    mod_L0_Mg <-
        MCMCglmm(Size ~ Sex + Days_Birth,
                 random = ~ animal + Year_fac + ID_Mom,
                 data = as.data.frame(tbl_L0),
                 pedigree = as.data.frame(tbl_ped |>
                                          select(id = ID_Juv,
                                                 dam = ID_Mom,
                                                 sire = ID_Dad)),
                 nitt = 5e5,
                 burnin = 1e4,
                 thin = 10)
    saveRDS(mod_L0_Mg, here("Output/animal_model_MCMCglmm_L0.rds"), compress = "xz")
    # Extremely similar to the brms output
}

# Loading the model from brms
mod_L0 <- readRDS(here("Output/animal_model_L0.rds"))

# Extracting the parameters of interest from the model
tbl_est_L0 <-
    bind_cols(
        get_sex_intercept(mod_L0),
        get_var_random(mod_L0)
    ) |>
    add_deriv_parameters()

# Difference between sexes
tbl_est_L0 |>
    select(Iter, contains("Int")) |>
    transmute(Diff = Int_SexM - Int_SexF) |>
    summarise(across(Diff, \(col_) summarise_chains(col_, with_p = TRUE)))

# Slope with birth date
mod_L0 |>
    fixef(summary = FALSE) |>
    as_tibble() |>
    summarise(across(Days_Birth, \(col_) summarise_chains(col_, with_p = TRUE)))


## ------------------------------------- Non-linear model for total growth ----

## Generating a dataset ready for the model (with multiple imputation)

format_nl <- function(df_) {
    df_ |>
        # Getting individuals with at least 2 records and known relatedness
        filter(n() >= 2, # NOTE Increasing to 3 does not improve CIs
               ID %in% rownames(mat_A),
               .by = ID) |>
        mutate(ID_PE = ID) |>  # Duplicating the ID column for permanent environment effect
        left_join(tbl_repro |>
                  select(ID = ID_Juv,
                         ID_Mom,
                         Date_Birth),
                  by = "ID") |>
        mutate(ID_Mom = if_else(is.na(ID_Mom), str_c("X_", ID), ID_Mom),
               .by = ID) |>
        mutate(Clutch_Size = map2_dbl(ID_Mom, Year,
                                      \(id, yr) { sum(tbl_repro[["ID_Mom"]] == id &
                                                      tbl_repro[["Year"]] == yr) }),
               Year       = year(Date_Birth),
               Year_fac   = as.factor(Year),
               Date_Birth = str_c(year(today()),
                                  month(Date_Birth),
                                  day(Date_Birth),
                                  sep = "-") |>
                            ymd(),
               Days_Birth = Date_Birth -
                            ymd(str_c(year(today()), "06", "01", sep = "-")))|>
        mutate(Days_Birth = Days_Birth - mean(Days_Birth),
               .by = Year)
}

# Formatting the datasets with and without imputation
lst_mult_nl <- map(lst_mult_capt, format_nl, .progress = TRUE)
tbl_nl_noimput <-
    format_nl(tbl_capt) |>
    filter(Size_NB_SE == 0)

# Subsetting the relatedness matrix
vec_ids_nl <- with(lst_mult_nl[[1]], unique(ID))
mat_A_nl <- mat_A[vec_ids_nl, vec_ids_nl]

## Running the non-linear brms model
# Prior distribution
prior_nl <-
    prior(normal(60, 10), nlpar = "Lmax") +
    prior(normal(0, 10), nlpar = "Lmax", coef = "SexM") +
    prior(uniform(0, 10), nlpar = "k") +
    prior(normal(0, 5), nlpar = "k", coef = "SexM")

# Setting up the formula
form_log <-
    bf(Size ~ Lmax  / (1 + ((Lmax - Size_NB) / Size_NB) * exp(- k * Age)),
       Lmax ~ Sex +
                  (1 | Year_fac) +
                  (0 + Sex || ID_PE) +
                  (0 + Sex || gr(ID, cov = A)) +
                  (1 | ID_Mom),
       k ~ Sex +
                  (1 | Year_fac) +
                  (0 + Sex || ID_PE) +
                  (0 + Sex || gr(ID, cov = A)) +
                  (1 | ID_Mom),
       nl = TRUE)

if (run_models_brms) {
    print("Running brms model for full growth")
    mod_nl <-
        brm_multiple(formula    = form_log,
                     data       = lst_mult_nl,
                     data2      = rep(list(list(A = as.matrix(mat_A_nl))),
                                      length(lst_mult_nl)),
                     chains     = n_chains_nl,
                     cores      = n_cores,
                     init       = rep(list(list(b_Lmax = array(data = c(68, 0)),
                                                b_k    = array(data = c(1, 0)))),
                                      n_chains_nl),
                     prior      = prior_nl,
                     seed       = seed,
                     save_pars  = save_pars(group = FALSE),
                     iter       = n_iter_nl,
                     thin       = thin_nl,
                     warmup     = n_wup_nl)
    saveRDS(mod_nl, here("Output/animal_model_nonlinear.rds"), compress = "xz")

    mod_nl_noimput <-
        brm(formula    = form_log,
            data       = tbl_nl_noimput,
            data2      = list(A = as.matrix(mat_A_nl)),
            chains     = n_chains_lm,
            cores      = n_cores,
            init       = rep(list(list(b_Lmax = array(data = c(68, 0)),
                                       b_k    = array(data = c(1, 0)))),
                             n_chains_lm),
            prior      = prior_nl,
            seed       = seed,
            save_pars  = save_pars(group = FALSE),
            iter       = n_iter_lm * 2,
            thin       = thin_lm * 2,
            warmup     = n_wup_lm)
    saveRDS(mod_nl_noimput, here("Output/animal_model_nonlinear_noimput.rds"), compress = "xz")
}

# Loading the model from brms
mod_NL <- readRDS(here("Output/animal_model_nonlinear.rds"))

# Extracting the parameters of interest from the model
tbl_est_NL <-
    full_join(
        get_sex_intercept(mod_NL, params = c("k", "Lmax")),
        get_var_random(mod_NL, params = c("k", "Lmax"))
    ) |>
    pivot_longer(contains("Sex"),
                 names_pattern = "(.*)_Sex([FM])",
                 names_to = c(".value", "Sex")) |>
    add_deriv_parameters(params = c("k", "Lmax")) |>
    select(Iter, Parameter, Sex, starts_with("Int"), starts_with("V"), everything())

# Difference in average between sexes
tbl_est_NL |>
    select(Iter, Parameter, Sex, Int) |>
    pivot_wider(names_from = Sex,
                names_prefix = "Int_",
                values_from = Int) |>
    mutate(Diff = Int_M - Int_F) |>
    summarise(across(Diff, \(col_) summarise_chains(col_, with_p = TRUE)),
              .by = Parameter) |>
    unpack(Diff)

# Log-ratio of V_A between sexes
tbl_est_NL |>
    select(Iter, Parameter, Sex, V_A) |>
    pivot_wider(names_from = Sex,
                names_prefix = "V_A_",
                values_from = V_A) |>
    mutate(LRV = log10(V_A_M/V_A_F)) |>
    summarise(across(LRV, \(col_) summarise_chains(col_, with_p = TRUE)),
              .by = Parameter) |>
    unpack(LRV)

# Log-ratio of h2 between sexes
tbl_est_NL |>
    select(Iter, Parameter, Sex, Herit) |>
    pivot_wider(names_from = Sex,
                names_prefix = "H_",
                values_from = Herit) |>
    mutate(LRV = log10(H_M/H_F)) |>
    summarise(across(LRV, \(col_) summarise_chains(col_, with_p = TRUE)),
              .by = Parameter) |>
    unpack(LRV)

## ------------------------------------------ Summary tables of the models ----

# Table of the output of mod_L0
smry_L0 <-
    tbl_est_L0 |>
    select(-starts_with("Int"), -Evolv) |>
    select(Iter, Mean, starts_with("V"), Herit) |>
    mutate(Parameter = "L0") |>
    pivot_wider(values_from = Mean:Herit,
                id_cols = "Iter",
                names_glue = "{Parameter}_{.value}",
                names_from = c(Parameter)) |>
    select(-Iter) |>
    map(summarise_chains) |>
    bind_rows(.id = "Parameter") |>
    separate(Parameter,
             into   = c("Parameter",  "Estimate"),
             sep    = "_",
             extra  = "merge")

write_ods(smry_L0,
          path  = here("Tables/Estimates_models.ods"),
          sheet = "L0")

# Table of the output of non-linear model
smry_nl <-
    tbl_est_NL |>
    select(-starts_with("Int"), -Evolv) |>
    select(Iter, Parameter, Sex, Mean, starts_with("V"), Herit) |>
    pivot_wider(values_from = Mean:Herit,
                id_cols = "Iter",
                names_glue = "{Parameter}_{Sex}_{.value}",
                names_from = c(Parameter, Sex)) |>
    select(-Iter) |>
    map(summarise_chains) |>
    bind_rows(.id = "Parameter") |>
    separate(Parameter,
             into   = c("Parameter", "Sex", "Estimate"),
             sep    = "_",
             extra  = "merge")

write_ods(smry_nl |>
          filter(Parameter == "k", Sex == "F") |>
          select(-Sex),
          path      = here("Tables/Estimates_models.ods"),
          append    = TRUE,
          sheet     = "k (F)")

write_ods(smry_nl |>
          filter(Parameter == "k", Sex == "M") |>
          select(-Sex),
          path      = here("Tables/Estimates_models.ods"),
          append    = TRUE,
          sheet     = "k (M)")

write_ods(smry_nl |>
          filter(Parameter == "Lmax", Sex == "F") |>
          select(-Sex),
          path      = here("Tables/Estimates_models.ods"),
          append    = TRUE,
          sheet     = "Lmax (F)")

write_ods(smry_nl |>
          filter(Parameter == "Lmax", Sex == "M") |>
          select(-Sex),
          path      = here("Tables/Estimates_models.ods"),
          append    = TRUE,
          sheet     = "Lmax (M)")

# Combined LaTeX tables
latex_smry <-
    bind_rows(smry_L0, smry_nl |> arrange(Parameter)) |>
    mutate(across(where(is.numeric),
                  \(col_) signif(col_, digits = 3)),
           across(where(is.numeric),
                  \(col_) {
                      if_else(abs(col_) < 1e-3,
                                  latex_sci(col_),
                                  as.character(col_))
                  }),
           Out = str_glue("${Median}$\\par $[{Low},{Up}]$")) |>
    select(Parameter, Sex, Estimate, Out) |>
    mutate(Estimate = factor(Estimate,
                             levels = c("Mean",
                                        "V_A",
                                        "V_M",
                                        "V_Y",
                                        "V_F",
                                        "V_PE",
                                        "V_R",
                                        "V_P",
                                        "Herit")),
           Parameter = factor(Parameter,
                              levels = c("L0", "k", "Lmax"))) |>
    arrange(Estimate, Parameter) |>
    mutate(Estimate = str_replace(Estimate,
                                  "V_([A-Z]+)",
                                  "$V_{\\\\text{\\1}}$") |>
                      str_replace("Mean",
                                  "$\\\\mu$") |>
                      str_replace("Herit",
                                  "$h^{2}$"),
           Parameter = recode(Parameter,
                              L0    = "$L_{0}$",
                              k     = "$k$",
                              Lmax  = "$L_{\\infty}$"),
           Sex = if_else(is.na(Sex), "B", Sex)) |>
    pivot_wider(names_from = c(Parameter, Sex),
                values_from = Out) |>
    rename_with(\(col_) str_replace(col_, "_(B)", "")) |>
    rename_with(\(col_) str_replace(col_, "_([FM])", " (\\1)")) |>
    mutate(across(everything(), \(col_) if_else(is.na(col_), "---", col_))) |>
    kable(format        = "latex",
          booktabs      = TRUE,
          digits        = 1,
          align         = "rllllll",
          col.names     = c("Estimate",
                            "Both",
                            "Female",
                            "Male",
                            "Female",
                            "Male"),
          linesep       = "",
          escape        = FALSE) |>
    add_header_above(header = c("",
                                "$L_{0}$" = 1,
                                "$k$" = 2,
                                "$L_{\\\\infty}$" = 2),
                     background = "blue!15!gray",
                     escape = FALSE,
                     line = FALSE,
                     line_sep = 0) |>
    kable_styling(latex_options = "striped",
                  stripe_color = "gray") |>
    column_spec(1,
                background = "blue!5!white!90!black",
                width = "1.2cm",
                latex_valign = "m") |>
    column_spec(c(2,3,4),
                width = "2.3cm",
                latex_valign = "m") |>
    row_spec(0, background = "blue!15!gray")

write_file(latex_smry, file = here("Tables/Estimates_models.tex"))

## ---------------------------------------- Graphical output of the models ----

## Gathering the output from all models
tbl_est_tot <-
    bind_rows(
        tbl_est_L0 |> mutate(Parameter = "L0", V_PE = V_R, V_R = NULL),
        tbl_est_NL |> select(-Iter)
    ) |>
    select(Parameter, Sex, Mean, starts_with("Int_"), starts_with("V_"), Herit, Evolv) |>
    mutate(Parameter = factor(Parameter, levels = c("L0", "k", "Lmax")),
           Sex = if_else(is.na(Sex), "B", Sex))

## Graphics for heritability of all parameters
p_h2 <-
    ggplot(tbl_est_tot) +
    geom_violin(aes(x = Parameter, y = Herit, fill = Sex),
                scale = "width") +
    stat_summary(aes(x = Parameter, y = Herit, group = paste(Parameter, Sex)),
                  fun = mean,
                  geom = "point",
                  position = position_dodge(width = 0.9),
                  size = 2) +
    xlab("Parameter") + ylab("Heritability")

cairo_pdf(here("Figures/Heritability.pdf"), width = 6, height = 6)
plot(p_h2)
dev.off()

## Graphics for evolvability of all parameters
p_I <-
    ggplot(tbl_est_tot) +
    geom_violin(aes(x = Parameter, y = sqrt(Evolv), fill = Sex),
                scale = "width") +
    stat_summary(aes(x = Parameter, y = sqrt(Evolv), group = paste(Parameter, Sex)),
                 fun = mean,
                 geom = "point",
                 position = position_dodge(width = 0.9),
                 size = 2) +
    xlab("Parameter") + ylab("Additive genetic coefficient of variation")

cairo_pdf(here("Figures/Evolvability.pdf"), width = 6, height = 6)
plot(p_I)
dev.off()

## Graphics for variance decomposition

# Setting a long table with all variances
tbl_var_long <-
    tbl_est_tot |>
    select(Parameter, Sex, starts_with("V_"), -V_P) |>
    pivot_longer(starts_with("V_"),
                 names_to = "Variance",
                 values_to = "Estimate") |>
    summarise(Estimate = median(Estimate),
              .by = c("Parameter", "Variance", "Sex")) |>
    mutate(Variance = factor(Variance, levels = c("V_F", "V_PE", "V_Y", "V_M", "V_A")))

p_var <-
    ggplot(tbl_var_long) +
    geom_col(aes(x = Parameter, y = Estimate, fill = Variance),
             position = "fill") +
    scale_fill_brewer(palette = "Set1", direction = -1) +
    facet_wrap(~ Sex, scale = "free")

cairo_pdf(here("Figures/Var_Decomp.pdf"), width = 6, height = 6)
plot(p_var)
dev.off()
