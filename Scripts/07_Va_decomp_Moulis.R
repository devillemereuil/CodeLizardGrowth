conflicted::conflicts_prefer(dplyr::filter)

##########################################################################################
##                        Quantitative genetics of the growth curves                    ##
##                                Analysis according to age                             ##
##                             "Genetics of growth analysis"                            ##
##                          Anaïs Aragon & Pierre de Villemereuil                       ##
##                                        (2023)                                        ##
##########################################################################################

## Packages required
library(tidyverse)
library(here)
library(brms)
library(furrr)
library(readODS)
library(patchwork)
library(kableExtra)
# remotes::install_github("devillemereuil/Reacnorm")
library(Reacnorm)

## Source the custom functions for this analysis
source("_functions.R")

## Setting up a (lucky) seed
seed <- 777
set.seed(seed)

## Setting up the ggplot2 theme
theme_set(theme_bw())
# Colours for males and females
col_sex <- c("#55aa00", "#00007f")
col_sex_points <- c("#a4ff7d", "#5560ff")

## Setting up parallelisation
n_cores  <- min(parallel::detectCores() - 2, 10)
options(mc.cores = n_cores)
plan(multicore)

## Setting up some parameters
# Number of samples to use for computing V_A / V_Gen per age
n_samples <- 1000
# Should the computation of the decomposition be run again?
run_decomp <- FALSE

## ----------------------------------- Loading the phenotypic dataset  ----

# Loading the dataset
tbl_capt <- 
    readRDS(here("Data/Moulis_capture_clean.rds")) |>
    mutate(Age = Age_Year)

## ----------------------------------- Loading the statistical models  ----

## Newborn model (L0)
# Loading the model
mod_L0 <- readRDS(here("Output/Moulis_animal_model_L0.rds"))
# Computing the estimates
tbl_est_L0 <-
    bind_cols(
        get_sex_intercept(mod_L0),
        get_var_random(mod_L0)
    ) |>
    add_deriv_parameters() |>
    mutate(Iter = 1:n(),
           Parameter = "L0",
           .before = 1) |>
    pivot_longer(contains("Int"),
                 names_to = "Sex",
                 values_to = "Int") |>
    mutate(Sex = str_remove(Sex, "Int_Sex"))

## Non-linear model
# Loading the model
mod_NL <- readRDS(here("Output/Moulis_animal_model_nonlinear.rds"))
# Computing the estimates
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
# Getting the residual variance of the non-linear model
vr_nl <-
    VarCorr(mod_NL, summary = FALSE)[["residual__"]][["sd"]]^2 |>
    as.vector()

## Grouping everything
tbl_est <-
    bind_rows(tbl_est_L0, tbl_est_NL) |>
    select(Iter, Parameter, Sex, Int, starts_with("V"), everything()) |>
    arrange(Iter)

## Subsetting the iterations
indices <- sample(1:max(tbl_est[["Iter"]]), size = n_samples)
tbl_est_sub <-
    tbl_est |>
    filter(Iter %in% indices)
vr_nl <- vr_nl[sort(indices)]

## -------- Generating the parameters to be used in the Reacnorm framework ----

# For this again, we will use the framework from de Villemereuil & Chevin (2024),
# implemented in the Reacnorm R package on GitHub (see above)

## Age vector
# The age values we will be computing the variances on
age_all <- 0:6  # All ages (almost)
age_rec <- 2:6  # Only recruits
# Computing the frequency of each age
wt_age_all  <-
    tbl_capt |>
    filter(Age >= 0, Age <= 6) |>
    count(Age) |>
    mutate(Prop = n / sum(n)) |>
    with(Prop)

wt_age_rec  <-
    tbl_capt |>
    filter(Age >= 2, Age <= 6) |>
    count(Age) |>
    mutate(Prop = n / sum(n)) |>
    with(Prop)

## Average L0 for Males and Females
size_L0 <-
    tbl_est_L0 |>
    select(Iter, Sex, Int) |>
    pivot_wider(names_from   = "Sex",
                values_from  = "Int") |>
    select(-Iter) |>
    colMeans()


## Computing the parameters intercept-value (by sex) and G-matrix
# Parameter intercepts for the females
theta_F <-
    tbl_est_sub |>
    filter(Sex == "F") |>
    select(Iter, Parameter, Int) |>
    pivot_wider(names_from = Parameter, values_from = Int) |>
    select(-Iter) |>
    apply(1, \(row_) { row_ }, simplify = FALSE)
names(theta_F) <- unique(tbl_est_sub[["Iter"]])

theta_M <-
    tbl_est_sub |>
    filter(Sex == "F") |>
    select(Iter, Parameter, Int) |>
    pivot_wider(names_from = Parameter, values_from = Int) |>
    select(-Iter) |>
    apply(1, \(row_) { row_ }, simplify = FALSE)
names(theta_M) <- unique(tbl_est_sub[["Iter"]])

# G-matrix for the parameters
G_F <-
    tbl_est_sub |>
    filter(Sex == "F") |>
    select(Iter, Parameter, V_A) |>
    pivot_wider(names_from = Parameter, values_from = V_A) |>
    select(-Iter) |>
    apply(1, \(row_) {
        G <- diag(row_);
        colnames(G) <- rownames(G) <- names(row_)
        return(G)
    }, simplify = FALSE)
names(G_F) <- unique(tbl_est_sub[["Iter"]])

G_M <-
    tbl_est_sub |>
    filter(Sex == "M") |>
    select(Iter, Parameter, V_A) |>
    pivot_wider(names_from = Parameter, values_from = V_A) |>
    select(-Iter) |>
    apply(1, \(row_) {
        G <- diag(row_);
        colnames(G) <- rownames(G) <- names(row_)
        return(G)
    }, simplify = FALSE)
names(G_M) <- unique(tbl_est_sub[["Iter"]])

# P-matrix for the parameters (for total variance
P_F <-
    tbl_est_sub |>
    filter(Sex == "F") |>
    select(Iter, Parameter, V_P) |>
    pivot_wider(names_from = Parameter, values_from = V_P) |>
    select(-Iter) |>
    apply(1, \(row_) {
        G <- diag(row_);
        colnames(G) <- rownames(G) <- names(row_)
        return(G)
    }, simplify = FALSE)
names(P_F) <- unique(tbl_est_sub[["Iter"]])

P_M <-
    tbl_est_sub |>
    filter(Sex == "M") |>
    select(Iter, Parameter, V_P) |>
    pivot_wider(names_from = Parameter, values_from = V_P) |>
    select(-Iter) |>
    apply(1, \(row_) {
        G <- diag(row_);
        colnames(G) <- rownames(G) <- names(row_)
        return(G)
    }, simplify = FALSE)
names(P_M) <- unique(tbl_est_sub[["Iter"]])

## Logistic growth functions and derivatives
# Female and Male growth curve functions
fn_logist <- expression(
    Lmax / (1 + ((Lmax - L0) / L0) * exp(-k * x))
)

## -------------------------- Computing V_Tot, V_Gen and V_A for each age  ----

## Compute V_tot for males and females
# Vectors of residual variance to set it at 0 for Age = 0
lst_vec_vr <-
    map(vr_nl, \(vr_) { c(0, rep(vr_, length(age_all) - 1)) })

if (run_decomp) {
    # Females
    var_tot_e_F <-
        future_pmap(list(theta_F, P_F, lst_vec_vr),
                    \(th, p, vr) { rn_vp_env(env      = age_all,
                                             shape    = fn_logist,
                                             theta    = th,
                                             V_theta  = p,
                                             var_res  = vr,
                                             width    = 20) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE)
    tbl_vtot_e_F <-
        var_tot_e_F |>
        map(\(vec) set_names(vec, age_all)) |>
        map(enframe, name = "Age", value = "V_Tot") |>
        list_rbind(names_to = "Iter") |>
        mutate(Sex = "F", Age = as.numeric(Age), .after = "Iter")
    
    # Males
    var_tot_e_M <-
        future_pmap(list(theta_M, P_M, lst_vec_vr),
                    \(th, p, vr) { rn_vp_env(env      = age_all,
                                             shape    = fn_logist,
                                             theta    = th,
                                             V_theta  = p,
                                             var_res  = vr,
                                             width    = 20) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE)
    tbl_vtot_e_M <-
        var_tot_e_M |>
        map(\(vec) set_names(vec, age_all)) |>
        map(enframe, name = "Age", value = "V_Tot") |>
        list_rbind(names_to = "Iter") |>
        mutate(Sex = "M", Age = as.numeric(Age), .after = "Iter")
    
    ## Compute V_Gen for males and females
    # Females
    var_gen_e_F <-
        future_map2(theta_F, G_F,
                    \(th, g) { rn_vgen(env      = age_all,
                                       shape    = fn_logist,
                                       theta    = th,
                                       G_theta  = g,
                                       average  = FALSE,
                                       width    = 20) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE)
    tbl_vgen_e_F <-
        var_gen_e_F |>
        map(\(vec) set_names(vec, age_all)) |>
        map(enframe, name = "Age", value = "V_Gen") |>
        list_rbind(names_to = "Iter") |>
        mutate(Sex = "F", Age = as.numeric(Age), .after = "Iter")
    
    # Males
    var_gen_e_M <-
        future_map2(theta_M, G_M,
                    \(th, g) { rn_vgen(env      = age_all,
                                       shape    = fn_logist,
                                       theta    = th,
                                       G_theta  = g,
                                       average  = FALSE,
                                       width    = 20) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE)
    tbl_vgen_e_M <-
        var_gen_e_M |>
        map(\(vec) set_names(vec, age_all)) |>
        map(enframe, name = "Age", value = "V_Gen") |>
        list_rbind(names_to = "Iter") |>
        mutate(Sex = "M", Age = as.numeric(Age), .after = "Iter")
    
    ## Compute V_A for males and females
    # Females
    var_a_e_nl_F <-
        future_map2(theta_F, G_F,
                    \(th, g) { rn_gamma_env(env      = age_all,
                                            shape    = fn_logist,
                                            theta    = th,
                                            G_theta  = g) |>
                              select(-matches(".*_.*_.*")) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE)
    tbl_va_e_F <-
        var_a_e_nl_F |>
        list_rbind(names_to = "Iter") |>
        as_tibble() |>
        select(!ends_with("Lmax_k")) |>
        rename(Age = Env) |>
        mutate(Sex = "F", .after = "Iter")
    
    # Males
    var_a_e_nl_M <-
        future_map2(theta_M, G_M,
                    \(th, g) { rn_gamma_env(env      = age_all,
                                            shape    = fn_logist,
                                            theta    = th,
                                            G_theta  = g) |>
                              select(-matches(".*_.*_.*")) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE)
    tbl_va_e_M <-
        var_a_e_nl_M |>
        list_rbind(names_to = "Iter") |>
        as_tibble() |>
        select(!ends_with("Lmax_k")) |>
        rename(Age = Env)|>
        mutate(Sex = "M", .after = "Iter")
    
    
    ## Concatening all the results
    tbl_vars <-
        bind_rows(tbl_vtot_e_F, tbl_vtot_e_M) |>
        full_join(bind_rows(tbl_vgen_e_F, tbl_vgen_e_M)) |>
        full_join(bind_rows(tbl_va_e_F, tbl_va_e_M))
    saveRDS(tbl_vars, here("Output/Moulis_Variances_Age_nonlinear.rds"))
} else {
    # Merely load the results
    tbl_vars <- readRDS(here("Output/Moulis_Variances_Age_nonlinear.rds"))
}

## --- Computing V_Add, V_A and V_AxE and their decomposition for all ages ----

## Getting average V_Tot for males and females to compute heritabilites
tbl_vtot_all <-
    tbl_vars |>
    summarise(V_Tot = weighted.mean(V_Tot[order(Age)], wt_age_all),
              .by = c(Iter, Sex)) |>
    pivot_wider(names_from = Sex, values_from = V_Tot) |>
    select(-Iter)

if (run_decomp) {
    ## Computing the full genetic decomposition for the females
    var_decomp_all_F <-
        future_map2(theta_F, G_F,
                    \(th, g) { rn_gen_decomp(env      = age_all,
                                             shape    = fn_logist,
                                             theta    = th,
                                             G_theta  = g,
                                             wt_env   = wt_age_all) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE) |>
        list_rbind() |>
        as_tibble() |>
        select(-matches(".*_.*_.*")) |>
        mutate(V_Tot = tbl_vtot_all[["F"]],
               Herit_Growth = V_Add / V_Tot,
               Herit_Size   = V_A / V_Tot,
               Herit_Shape    = V_AxE / V_Tot)
    
    saveRDS(var_decomp_all_F, file = here("Output/Moulis_var_decomp_all_F.rds"))

    ## Computing the full genetic decomposition for the males
    var_decomp_all_M <-
        future_map2(theta_M, G_M,
                    \(th, g) { rn_gen_decomp(env      = age_all,
                                             shape    = fn_logist,
                                             theta    = th,
                                             G_theta  = g,
                                             wt_env   = wt_age_all) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE) |>
        list_rbind() |>
        as_tibble() |>
        select(-matches(".*_.*_.*")) |>
        mutate(V_Tot = tbl_vtot_all[["M"]],
               Herit_Growth = V_Add / V_Tot,
               Herit_Size   = V_A / V_Tot,
               Herit_Shape    = V_AxE / V_Tot)
    
    saveRDS(var_decomp_all_M, file = here("Output/Moulis_var_decomp_all_M.rds"))

    write_ods(var_decomp_all_F |>
              pivot_longer(everything(),
                           names_to = "Parameter",
                           values_to = "Value") |>
              summarise(Mean = mean(Value) |> signif(digits = 3),
                        Median = median(Value) |> signif(digits = 3),
                        Low = quantile(Value, probs = 0.025)[1]  |> signif(digits = 3),
                        Up = quantile(Value, probs = 0.975)[1] |> signif(digits = 3),
                        .by = Parameter),
              path = here("Tables/Moulis_VA_decomp.ods"),
              sheet = "Females (all)")
    write_ods(var_decomp_all_M |>
              pivot_longer(everything(),
                           names_to = "Parameter",
                           values_to = "Value") |>
              summarise(Mean = mean(Value) |> signif(digits = 3),
                        Median = median(Value) |> signif(digits = 3),
                        Low = quantile(Value, probs = 0.025)[1]  |> signif(digits = 3),
                        Up = quantile(Value, probs = 0.975)[1] |> signif(digits = 3),
                        .by = Parameter),
              path = here("Tables/Moulis_VA_decomp.ods"),
              append = TRUE,
              sheet = "Males (all)")
} else {
    # Merely load the results
    var_decomp_all_F <- readRDS(here("Output/Moulis_var_decomp_all_F.rds"))
    var_decomp_all_M <- readRDS(here("Output/Moulis_var_decomp_all_M.rds"))
}

# Testing sex differences
tbl_contrast <-
    bind_rows(var_decomp_all_F |> mutate(Sex = "F", .before = 1),
          var_decomp_all_M |> mutate(Sex = "M", .before = 1)) |>
    reframe(across(-Sex, \(vec_) {log10(vec_[Sex == "M"]/vec_[Sex == "F"])})) |>
    map(\(col_) summarise_chains(col_, with_p = TRUE)) |>
    bind_rows(.id = "Param")

## --- Computing V_Add, V_A and V_AxE and their decomposition for recruits ----

## Getting average V_Tot for males and females to compute heritabilites
tbl_vtot_rec <-
    tbl_vars |>
    filter(Age %in% age_rec) |>
    summarise(V_Tot = weighted.mean(V_Tot[order(Age)], wt_age_rec),
              .by = c(Iter, Sex)) |>
    pivot_wider(names_from = Sex, values_from = V_Tot) |>
    select(-Iter)

if (run_decomp) {
    ## Computing the full genetic decomposition for the females
    var_decomp_rec_F <-
        future_map2(theta_F, G_F,
                    \(th, g) { rn_gen_decomp(env      = age_rec,
                                             shape    = fn_logist,
                                             theta    = th,
                                             G_theta  = g,
                                             wt_env   = wt_age_rec) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE) |>
        list_rbind() |>
        as_tibble() |>
        select(-matches(".*_.*_.*")) |>
        mutate(V_Tot = tbl_vtot_rec[["F"]],
               Herit_Growth = V_Add / V_Tot,
               Herit_Size   = V_A / V_Tot,
               Herit_Shape    = V_AxE / V_Tot)

    saveRDS(var_decomp_rec_F, file = here("Output/Moulis_var_decomp_rec_F.rds"))

    ## Computing the full genetic decomposition for the males
    var_decomp_rec_M <-
        future_map2(theta_M, G_M,
                    \(th, g) { rn_gen_decomp(env      = age_rec,
                                             shape    = fn_logist,
                                             theta    = th,
                                             G_theta  = g,
                                             wt_env   = wt_age_rec) },
                    .options = furrr_options(seed = TRUE),
                    .progress = TRUE) |>
        list_rbind() |>
        as_tibble() |>
        select(-matches(".*_.*_.*")) |>
        mutate(V_Tot = tbl_vtot_rec[["M"]],
               Herit_Growth = V_Add / V_Tot,
               Herit_Size   = V_A / V_Tot,
               Herit_Shape    = V_AxE / V_Tot)

    saveRDS(var_decomp_rec_M, file = here("Output/Moulis_var_decomp_rec_M.rds"))

    write_ods(var_decomp_rec_F |>
              pivot_longer(everything(),
                           names_to = "Parameter",
                           values_to = "Value") |>
              summarise(Mean = mean(Value) |> signif(digits = 3),
                        Median = median(Value) |> signif(digits = 3),
                        Low = quantile(Value, probs = 0.025)[1]  |> signif(digits = 3),
                        Up = quantile(Value, probs = 0.975)[1] |> signif(digits = 3),
                        .by = Parameter),
              path = here("Tables/Moulis_VA_decomp.ods"),
              append = TRUE,
              sheet = "Females (recruits)")
    write_ods(var_decomp_rec_M |>
              pivot_longer(everything(),
                           names_to = "Parameter",
                           values_to = "Value") |>
              summarise(Mean = mean(Value) |> signif(digits = 3),
                        Median = median(Value) |> signif(digits = 3),
                        Low = quantile(Value, probs = 0.025)[1]  |> signif(digits = 3),
                        Up = quantile(Value, probs = 0.975)[1] |> signif(digits = 3),
                        .by = Parameter),
              path = here("Tables/Moulis_VA_decomp.ods"),
              append = TRUE,
              sheet = "Males (recruits)")
} else {
    # Merely load the results
    var_decomp_rec_F <- readRDS(here("Output/Moulis_var_decomp_rec_F.rds"))
    var_decomp_rec_M <- readRDS(here("Output/Moulis_var_decomp_rec_M.rds"))
}

## ------------------------------------------ Summary table for estimates  ----

# Format the estimates
tbl_est <-
    list(All_F      = var_decomp_all_F,
         All_M      = var_decomp_all_M,
         Recruit_F  = var_decomp_rec_F,
         Recruit_M  = var_decomp_rec_M) |>
    bind_rows(.id = "Context") |>
    pivot_longer(-Context,
                 names_to = "Parameter",
                 values_to = "Value") |>
    summarise(Mean = mean(Value) |> signif(digits = 3),
              Median = median(Value) |> signif(digits = 3),
              Low = quantile(Value, probs = 0.025)[1]  |> signif(digits = 3),
              Up = quantile(Value, probs = 0.975)[1] |> signif(digits = 3),
              .by = c(Parameter, Context)) |>
    separate(Context,
             into = c("Stages", "Sex")) |> 
    mutate(Parameter = str_replace(Parameter, "V_AxE", "V_AxAge") |>
                       str_replace("Herit_SxA", "Herit_Shape"),
           Estimate = factor(Parameter,
                              levels = c("V_Add", "V_A", "V_AxAge", "V_Tot",
                                         "Herit_Growth", "Herit_Size", "Herit_Shape",
                                         "Gamma_L0", "Gamma_k", "Gamma_Lmax",
                                         "Iota_L0", "Iota_k", "Iota_Lmax"))) |>
    arrange(Stages, Sex, Estimate) |>
    mutate(across(where(is.numeric),
                  \(col_) signif(col_, digits = 3)),
           across(where(is.numeric),
                  \(col_) {
                      if_else(col_ < 1e-3,
                                  latex_sci(col_),
                                  as.character(col_))
                  }),
           Out = str_glue("${Median}$\\par $[{Low},{Up}]$")) |>
    select(Stages, Sex, Estimate, Out) |>
    pivot_wider(names_from = c(Stages, Sex),
                values_from = Out)

latex_smry <-
    tbl_est |>
    mutate(Estimate = str_replace(Estimate,
                                  "V_([a-zA-Z]+)",
                                  "$V_{\\\\text{\\1}}$") |>
                      str_replace("Herit_([a-zA-Z]+)",
                                  "$h^2_{\\\\text{\\1}}$") |>
                      str_replace("Gamma_(.*)", "$\\\\gamma_{\\1}$") |>
                      str_replace("Iota_(.*)", "$\\\\iota_{\\1}$") |>
                      str_replace("L0", "L_{0}") |>
                      str_replace("Lmax", "L_{\\\\infty}")) |>
    kable(format        = "latex",
          booktabs      = TRUE,
          digits        = 1,
          align         = "rllll",
          linesep       = "",
          escape        = FALSE,
          col.names     = case_when(str_detect(colnames(tbl_est), "_F") ~ "Females",
                                    str_detect(colnames(tbl_est), "_M") ~ "Males",
                                    .default = colnames(tbl_est))) |>
    kable_styling(latex_options = "striped",
                  stripe_color = "gray") |>
    column_spec(1,
                background = "blue!5!white!90!black",
                width = "1.2cm",
                latex_valign = "m") |>
    column_spec(2:5,
                width = "2.3cm",
                latex_valign = "m") |>
    row_spec(0, background = "blue!15!gray") |>
    add_header_above(c("", "All ages" = 2, "Adults only" = 2),
                     background = "blue!15!gray",
                     line = FALSE,
                     line_sep = 0) |>
    pack_rows("Variances", 1, 4,
              background = "blue!5!white!90!black",
              escape = FALSE) |>
    pack_rows("Heritabilities", 5, 7,
              background = "blue!5!white!90!black",
              escape = FALSE) |>
    pack_rows("$\\gamma$-decomposition", 8, 10,
              background = "blue!5!white!90!black",
              escape = FALSE) |>
    pack_rows("$\\iota$-decomposition", 11, 13,
              background = "blue!5!white!90!black",
              escape = FALSE)

write_file(latex_smry, file = here("Tables/Moulis_Estimates_reacnorm.tex"))

## ------------------------------------------ Graphics output against age  ----

## Variances against age
tbl_plot_vars <-
    tbl_vars |>
    select(!starts_with("Gamma"), -V_Tot) |>
    mutate(Sex = if_else(Age == 0, "J", Sex)) |>
    pivot_longer(starts_with("V_"),
                 names_to = "Variance",
                 values_to = "Estimate") |>
    mutate(Age_fac = factor(Age, levels = age_all),
           Variance = recode(Variance, V_A = "V[A]", V_Gen = "V[Gen]"),
           Sex = factor(Sex, levels = c("F", "M", "J")))

p_vars <-
    ggplot(tbl_plot_vars) +
    geom_violin(aes(x = Age_fac, y = Estimate, fill = Sex),
                scale = "width") +
    stat_summary(aes(x = Age_fac, y = Estimate, colour = Sex),
                 geom = "point",
                 fun = "mean",
                 position = position_dodge(width = 1),
                 size = 2) +
    facet_wrap(~ Variance, labeller = label_parsed) +
    scale_colour_manual(values = c(col_sex, "black")) +
    scale_fill_manual(values = c(col_sex_points, "#bababa")) +
    xlab("Age") + ylab("Variance")

cairo_pdf(here("Figures/Moulis_Genetic_variance_Age.pdf"), width = 8, height = 6)
plot(p_vars)
dev.off()

## Heritability against age
tbl_plot_h2 <-
    tbl_vars |>
    transmute(Sex       = Sex,
              Age       = Age,
              Herit     = V_A / V_Tot,
              Herit_B   = V_Gen / V_Tot) |>
    mutate(Sex = if_else(Age == 0, "J", Sex)) |>
    pivot_longer(starts_with("Herit"),
                 names_to = "Herit",
                 values_to = "Estimate") |>
    mutate(Age_fac = factor(Age, levels = age_all),
           Herit    = recode(Herit,
                             Herit = "h^2",
                             Herit_B = "H^2"),
           Sex      = factor(Sex, levels = c("F", "M", "J")))

p_herit <-
    ggplot(tbl_plot_h2) +
    geom_violin(aes(x = Age_fac, y = Estimate, fill = Sex),
                scale = "width") +
    stat_summary(aes(x = Age_fac, y = Estimate, colour = Sex),
                 geom = "point",
                 fun = "mean",
                 position = position_dodge(width = 1),
                 size = 2) +
    facet_wrap(~ Herit, labeller = label_parsed) +
    scale_colour_manual(values = c(col_sex, "black")) +
    scale_fill_manual(values = c(col_sex_points, "#bababa")) +
    xlab("Age") + ylab("Heritability")

cairo_pdf(here("Figures/Moulis_Herit_Age.pdf"), width = 8, height = 6)
plot(p_herit)
dev.off()


## Both V_A / h² against age
tbl_plot_both <-
    bind_rows(tbl_plot_vars |> filter(Variance == "V[A]") |> rename(Parameter = Variance),
              tbl_plot_h2 |> filter(Herit == "h^2") |> rename(Parameter = Herit)) |>
    mutate(Parameter = str_c("italic(", Parameter, ")") |>
                       as_factor(),
           Sex = recode(Sex, M = "Male", F = "Female", J = "Juvenile"))


p_both <-
    ggplot(tbl_plot_both) +
    geom_violin(data    = tbl_plot_both |>
                          filter(between(Estimate,
                                         quantile(Estimate, probs = 0.025),
                                         quantile(Estimate, probs = 0.975)),
                                 .by = c(Age_fac, Parameter, Sex)),
                mapping = aes(x = Age_fac, y = Estimate, fill = Sex, colour = Sex),
                linewidth = 0,
                scale   = "width") +
    geom_violin(aes(x = Age_fac, y = Estimate, colour = Sex),
                fill = "#FFFFFF00",
                scale = "width") +
    stat_summary(aes(x = Age_fac, y = Estimate, colour = Sex),
                 geom = "point",
                 fun = "mean",
                 position = position_dodge(width = 1),
                 size = 2) +
    facet_wrap(~ Parameter, ncol =1, labeller = label_parsed, scale = "free_y") +
    scale_colour_manual(name = "Sex/Stage",
                        values = c(col_sex, "black")) +
    scale_fill_manual(name = "Sex/Stage",
                      values = c(col_sex_points, "#bababa")) +
    xlab("Age") + ylab("Value")

cairo_pdf(here("Figures/Moulis_V_A_h2_Age.pdf"), width = 6, height = 8)
plot(p_both)
dev.off()


## Gamma against age
tbl_plot_g <-
    tbl_vars |>
    select(!starts_with("V_")) |>
    pivot_longer(starts_with("Gamma_"),
                 names_to = "Gamma",
                 values_to = "Estimate") |>
    mutate(Age_fac = factor(Age, levels = age_all),
           Gamma = Gamma |>
                   str_remove("Gamma_") |>
                   recode(lmax = "L[infinity]"),
           Sex = factor(Sex, levels = c("F", "M"))) |>
    summarise(Low = quantile(Estimate, probs = 0.05),
              Up = quantile(Estimate, probs = 0.95),
              Middle = mean(Estimate),
              .by = c(Sex, Age, Gamma, Age_fac)) |>
    mutate(Gamma = factor(Gamma,
                          levels = c("L0", "k", "Lmax")) |>
                   recode(L0    = "italic(L[0])",
                          Lmax  = "italic(L[infinity])",
                          k     = "italic(k)"),
           Sex   = recode(Sex, M = "Male", F = "Female"),
           Title = "γ-decomposition")

p_gamma <-
    ggplot(tbl_plot_g) +
    geom_ribbon(aes(x = Age, ymin = Low, ymax = Up, fill = Sex, group = paste0(Gamma, Sex)),
                alpha = 0.2) +
    geom_line(aes(x = Age, y = Middle, colour = Sex, linetype = Gamma)) +
    facet_wrap(~ Title) +
    scale_colour_manual(values = col_sex) +
    scale_fill_manual(values = col_sex_points) +
    scale_linetype_manual(name = "Parameter",
                          values = c("dotted", "solid", "longdash"),
                          labels = scales::label_parse()) +
    xlab("Age") + ylab("")

cairo_pdf(here("Figures/Moulis_Gamma_Age.pdf"), width = 7, height = 6)
plot(p_gamma)
dev.off()

## Full graph for the decomposition

cairo_pdf(here("Figures/Moulis_Full_decomp.pdf"), width = 10, height = 6)
plot(p_both + p_gamma + plot_layout(guides = "collect") & theme(legend.position = "top"))
dev.off()
