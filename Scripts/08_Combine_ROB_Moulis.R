conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(brms::ar)
conflicted::conflicts_prefer(dplyr::lag)
conflicted::conflicts_prefer(kableExtra::group_rows)

##########################################################################################
##                     Study and selection of the different growth curves               ##
##                                   (Moulis version)                                   ##
##                             "Genetics of growth analysis"                            ##
##                          Anaïs Aragon & Pierre de Villemereuil                       ##
##                                        (2024)                                        ##
##########################################################################################

## Packages required
library(tidyverse)
library(here)
library(brms)
library(future)
library(readODS)
library(patchwork)
library(kableExtra)
library(scales)

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

## ----------------------------------- Growth curve selection (WAIC table) ----

# Loading the separate WAIC tables
waic_rob <- read_ods(here("Tables/WAIC.ods"))
waic_mls <- read_ods(here("Tables/Moulis_WAIC.ods"))

# Combining the tables
both_waic <-
    both_waic <-
    left_join(waic_rob, waic_mls,
              by     = "Model",
              suffix = c("_ROB", "_Moulis"))

# Combined summary table
tex_waic <-
    both_waic |>
    kable(format        = "latex",
          booktabs      = TRUE,
          digits        = 1,
          escape        = FALSE,
          col.names     = colnames(both_waic) |>
                          str_replace("Delta_WAIC", "$\\\\Delta$WAIC") |>
                          str_replace("SE", "Std. Error") |>
                          str_remove("_(ROB|Moulis)")) |>
    column_spec(1, bold = TRUE, background = "blue!5!white!90!black") |>
    row_spec(0, bold = TRUE, background = "blue!15!gray") |>
    add_header_above(c("", "Wild" = 3, "Mesocosm" = 3),
                     background = "blue!15!gray",
                     bold = TRUE,
                     line = FALSE,
                     line_sep = 0) |>
    kable_styling() |>
    str_remove_all("\\\\(begin|end)\\{table\\}(\n)?")

write_file(tex_waic, file = here("Tables/WAIC_both.tex"))

## --------------------------------------------------- Final growth curves ----

# Loading the graphs
p_rob <- readRDS(here("Output/Object_Final_Growth_Curve_ROB.rds"))
p_mls <- readRDS(here("Output/Object_Final_Growth_Curve_Moulis.rds"))

# Combined graphics
cairo_pdf(here("Figures/Both_Final_Growth_Curve.pdf"), width = 12, height = 5)
plot((p_rob + ggtitle("Wild")) + (p_mls + ggtitle("Mesocosm")) +
     plot_layout(guides = "collect") &
     theme(legend.position = "bottom",
           legend.spacing.x = unit(3, "cm")))
dev.off()

## ------------------------------------------------------- Model estimates ----

# List of parameters
vec_pars <- c("L0", "k (F)", "k (M)", "Lmax (F)", "Lmax (M)")
names(vec_pars) <- vec_pars

# Loading estimates from ROB
tbl_est_rob <-
    map(vec_pars,
        \(p_) {
            read_ods(here("Tables/Estimates_models.ods"),
                     sheet = p_)
        }) |>
    bind_rows(.id = "Sheet") |>
    mutate(Population = "ROB")

# Loading estimates from Moulis
tbl_est_mls <-
    map(vec_pars,
        \(p_) {
          read_ods(here("Tables/Moulis_Estimates_models.ods"),
                   sheet = p_)
        }) |>
    bind_rows(.id = "Sheet") |>
    mutate(Population = "Moulis")

# Combining and formatting the tbl
tbl_est <-
    bind_rows(tbl_est_rob, tbl_est_mls) |>
    mutate(Sex = case_when(str_detect(Sheet, "\\(F\\)")  ~ "F",
                           str_detect(Sheet, "\\(M\\)")  ~ "M",
                           .default = "B"),
           Sheet = NULL,
           Population = as_factor(Population),
           across(where(is.numeric),
                  \(col_) signif(col_, digits = 3)),
           across(where(is.numeric),
                  \(col_) {
                      if_else(col_ < 1e-2,
                                  latex_sci(col_, 3),
                                  to_signif(col_, 3))
                  }),
           Out = str_glue("${Median}$\\par\\scriptsize $[{Low},{Up}]$")) |>  
    select(Population, Sex, Parameter, Estimate, Out) |>
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
    arrange(Population, Estimate, Parameter) |>
    mutate(Estimate = str_replace(Estimate,
                                  "V_([A-Z]+)",
                                  "$V_{\\\\text{\\1}}$") |>
                      str_replace("Mean",
                                  "$\\\\mu$") |>
                      str_replace("Herit",
                                  "$h^{2}$") |>
                      str_replace("V_\\{\\\\text\\{P\\}\\}",
                                  "V_{\\\\text{Par}}"),
           Parameter = recode(Parameter,
                              L0    = "$L_{0}$",
                              k     = "$k$",
                              Lmax  = "$L_{\\\\infty}$")) |>
    pivot_wider(names_from = c(Population, Parameter, Sex),
                values_from = Out) |>
    mutate(across(everything(), \(col_) if_else(is.na(col_), "---", col_)))

# Combined LaTeX tables
latex_est <-
    tbl_est |>
    kable(format        = "latex",
          booktabs      = TRUE,
          digits        = 1,
          align         = "rcccccccccc",
          linesep       = "",
          escape        = FALSE,
          col.names     = case_when(
              colnames(tbl_est) == "Estimate"         ~ "Estim.",
              str_detect(colnames(tbl_est), "_B$")    ~ "",
              str_detect(colnames(tbl_est), "_F$")    ~ "Female",
              str_detect(colnames(tbl_est), "_M$")    ~ "Male"
          )) |>
    kable_styling(latex_table_env = NULL,
                  latex_options = "striped",
                  stripe_color = "gray") |>
    column_spec(1,
                background = "blue!5!white!90!black",
                width = "0.8cm",
                latex_valign = "m") |>
    column_spec(2:11,
                width = "1.7cm",
                latex_valign = "m") |>
    row_spec(0, background = "blue!15!gray") |>
    add_header_above(colnames(tbl_est) |>
                     str_remove("(ROB|Moulis)_") |>
                     str_remove("Estimate") |>
                     str_remove("_[BFM]") |>
                     group_col(),
                     background = "blue!15!gray",
                     escape = FALSE,
                     line = FALSE,
                     line_sep = 0) |>
    add_header_above(case_when(str_detect(colnames(tbl_est), "ROB")    ~ "Wild",
                               str_detect(colnames(tbl_est), "Moulis") ~ "Mesocosm",
                               .default = "") |>
                     group_col(),
                     background = "blue!15!gray",
                     bold = TRUE,
                     line = TRUE,
                     line_sep = 3) |>
    str_remove_all("\\\\(begin|end)\\{table\\}(\n)?") |>
    str_replace("\\\\cmidrule",
                "\\\\arrayrulecolor{blue!15!gray} \\\\specialrule{6pt}{0pt}{-6pt} \\\\arrayrulecolor{black} \\\\cmidrule") |>
    str_replace(fixed("\\cellcolor{blue!15!gray}{Estim.}"),
                "\\arrayrulecolor{blue!15!gray} \\specialrule{6pt}{0pt}{-6pt} \\arrayrulecolor{black} \\cmidrule(l{3pt}r{3pt}){3-4} \\cmidrule(l{3pt}r{3pt}){5-6} \\cmidrule(l{3pt}r{3pt}){8-9} \\cmidrule(l{3pt}r{3pt}){10-11}\n \\cellcolor{blue!15!gray}{Estim.}")

write_file(latex_est, file = here("Tables/Both_estimates_models.tex"))

## ----------------------------------------------------- V_A Decomposition ----

# List of parameters
vec_context <- c(F = "Females (all)", M = "Males (all)")

# Loading decomposition from ROB
tbl_decomp_rob <-
    map(vec_context,
        \(p_) {
          read_ods(here("Tables/VA_decomp.ods"),
                   sheet = p_)
        }) |>
    bind_rows(.id = "Sex") |>
    mutate(Population = "ROB")

# Loading decomposition from Moulis
tbl_decomp_mls <-
    map(vec_context,
        \(p_) {
          read_ods(here("Tables/Moulis_VA_decomp.ods"),
                   sheet = p_)
        }) |>
    bind_rows(.id = "Sex") |>
    mutate(Population = "Moulis")

# Combining and formatting the tbl
tbl_decomp <-
    bind_rows(tbl_decomp_rob, tbl_decomp_mls) |>
    mutate(Population = as_factor(Population),
           Parameter = str_replace(Parameter, "V_AxE", "V_AxAge") |>
                       str_replace("Herit_SxA", "Herit_Shape"),
           Estimate = factor(Parameter,
                              levels = c("V_Add", "V_A", "V_AxAge", "V_Tot",
                                         "Herit_Growth", "Herit_Size", "Herit_Shape",
                                         "Gamma_L0", "Gamma_k", "Gamma_Lmax",
                                         "Iota_L0", "Iota_k", "Iota_Lmax"))) |>
    arrange(Population, Sex, Estimate) |>
    mutate(across(where(is.numeric),
                  \(col_) signif(col_, digits = 3)),
           across(where(is.numeric),
                  \(col_) {
                      if_else(col_ < 1e-3,
                                  latex_sci(col_, 3),
                                  to_signif(col_, 3))
                  }),
           Out = str_glue("${Median}$\\par\\scriptsize $[{Low},{Up}]$")) |>
    select(Population, Sex, Estimate, Out) |>
    pivot_wider(names_from = c(Population, Sex),
                values_from = Out)

# Combined LaTeX tables
latex_decomp <-
    tbl_decomp |>
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
          align         = "rcccc",
          linesep       = "",
          escape        = FALSE,
          col.names     = colnames(tbl_decomp) |>
                          str_remove("(ROB|Moulis)_") |>
                          str_replace("^F$", "Females") |>
                          str_replace("^M$", "Males") |>
                          str_replace("Estimate", "Estim.")) |>
    kable_styling(latex_table_env = NULL,
                  latex_options = "striped",
                  stripe_color = "gray") |>
    column_spec(1,
                background = "blue!5!white!90!black",
                width = "1cm",
                latex_valign = "m") |>
    column_spec(2:5,
                width = "2.3cm",
                latex_valign = "m") |>
    row_spec(0, background = "blue!15!gray") |>
    add_header_above(c("", "Wild" = 2, "Mesocosm" = 2),
                     background = "blue!15!gray",
                     bold = TRUE,
                     line = TRUE,
                     line_sep = 3) |>
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
              escape = FALSE) |>
    str_remove_all("\\\\(begin|end)\\{table\\}(\n)?") |>
    str_replace("\\\\cmidrule",
                "\\\\arrayrulecolor{blue!15!gray} \\\\specialrule{6pt}{0pt}{-6pt} \\\\arrayrulecolor{black} \\\\cmidrule")

write_file(latex_decomp, file = here("Tables/Both_estimates_reacnorm.tex"))

## ------------------------------------------------------------ V_A by Age ----

# Loading the variance decomposition against age
tbl_decomp_rob <-
    readRDS(here("Output/Variances_Age_nonlinear.rds")) |>
    mutate(Population = "ROB")
tbl_decomp_mls <-
    readRDS(here("Output/Moulis_Variances_Age_nonlinear.rds")) |>
    mutate(Population = "Moulis")

# Combining the estimates against age
tbl_decomp <-
    bind_rows(tbl_decomp_rob, tbl_decomp_mls) |>
    mutate(Population = as_factor(Population),
           Herit = V_A / V_Tot) |>
    select(-V_Gen, -V_Tot) |>
    pivot_longer(c(starts_with("Gamma_"), starts_with("V_"), "Herit"),
                 names_to = "Parameter",
                 values_to = "Value") |>
    mutate(Type = case_when(Parameter == "V_A"              ~ "V_A",
                            Parameter == "Herit"            ~ "Herit",
                            str_detect(Parameter, "Gamma")  ~ "Gamma") |>
                  factor(levels = c("V_A", "Herit", "Gamma")),
           Parameter = factor(Parameter,
                              levels = c("V_A", "Herit",
                                         "Gamma_L0", "Gamma_k", "Gamma_Lmax"))) |>
    arrange(Population, Sex, Age, Type, Parameter) |>
    mutate(Age_fac = factor(Age, levels = 0:6),
           Type = recode(Type,
                         V_A    = "italic(V[A])",
                         Herit  = "italic(h^2)",
                         Gamma  = "italic(γ)-decomposition"),
           Parameter = recode(Parameter,
                              V_A        = "italic(V[A])",
                              Herit      = "italic(h^2)",
                              Gamma_L0   = "italic(γ[L[0]])",
                              Gamma_k    = "italic(γ[k])",
                              Gamma_Lmax = "italic(γ[L[infinity]])"),
           Sex = recode(Sex, M = "Males", F = "Females"))

generate_plot_decomp <- function(tbl_, title, ylab) {
    tbl_nog_ <-
        tbl_ |>
        filter(!str_detect(Type, "γ"))

    tbl_g_ <-
        tbl_ |>
        filter(str_detect(Type, "γ")) |>
        summarise(Middle = median(Value),
                  Low    = quantile(Value, probs = 0.05),
                  Up     = quantile(Value, probs = 0.95),
                  .by = c(Age, Age_fac, Sex, Parameter, Type))

    p <-
        ggplot(tbl_) +
        geom_violin(data = tbl_nog_ |>
                           filter(between(Value,
                                          quantile(Value, probs = 0.05),
                                          quantile(Value, probs = 0.95)),
                                  .by = c(Age_fac, Parameter, Sex, Type)),
                    mapping = aes(x         = Age_fac,
                                  y         = Value,
                                  fill      = Sex,
                                  colour    = Sex),
                    linewidth = 0,
                    scale   = "width") +
        geom_violin(data = tbl_nog_,
                    mapping = aes(x         = Age_fac,
                                  y         = Value,
                                  group     = paste0(Age_fac, Sex)),
                    fill    = "#FFFFFF00",
                    scale   = "width") +
        stat_summary(data = tbl_nog_,
                     mapping = aes(x        = Age_fac,
                                   y        = Value,
                                   colour   = Sex),
                     geom = "point",
                     fun = "mean",
                     position = position_dodge(width = 1),
                     size = 2) +
        geom_ribbon(data = tbl_g_,
                    mapping = aes(x         = Age + 1,
                                  ymin      = Low,
                                  ymax      = Up,
                                  fill      = Sex,
                                  group     = paste0(Parameter, Sex)),
                    alpha = 0.2) +
        geom_line(data = tbl_g_,
                  mapping = aes(x           = Age + 1,
                                y           = Middle,
                                colour      = Sex,
                                linetype    = Parameter)) +
        facet_grid(Type ~ .,
                   labeller = label_parsed,
                   scale = "free") +
        scale_colour_manual(name = "Sex",
                            values = c(col_sex, "black")) +
        scale_fill_manual(name = "Sex",
                          values = c(col_sex_points, "#bababa")) +
        scale_linetype_manual(name = "Parameter",
                              values = c("dotted", "solid", "longdash"),
                              labels = scales::label_parse()) +
         xlab("Age") + ylab(ylab) +
         ggtitle(title)

    return(p)
}

p_decomp_rob <-
    generate_plot_decomp(tbl_decomp |> filter(Population == "ROB"),
                     title = "Wild",
                     ylab  = "Value")
p_decomp_mls <-
    generate_plot_decomp(tbl_decomp |> filter(Population == "Moulis"),
                         title = "Mesocosm",
                         ylab  = "")

cairo_pdf(here("Figures/Both_full_decomp.pdf"), width = 10, height = 10)
plot(p_decomp_rob +
     p_decomp_mls +
     plot_layout(guides = "collect", widths = c(1, 1)) &
     theme(legend.position = "bottom"))
dev.off()

