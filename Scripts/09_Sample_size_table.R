conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::lag)

##########################################################################################
##                             Generating a sample size table                           ##
##                             "Genetics of growth analysis"                            ##
##                          Anaïs Aragon & Pierre de Villemereuil                       ##
##                                        (2026)                                        ##
##########################################################################################

## Packages required
library(tidyverse)
library(here)
library(brms)
library(kableExtra)

## Source the custom functions for this analysis
source("_functions.R")

## Setting up a (lucky) seed
seed <- 777
set.seed(seed)

## -------------------------------------------------- Setting up the table ----

# This table will contain the levels envisioned for the Sample Size Table
# in the article
tbl_samp <- crossing(Population = c("ROB", "Moulis"),
                     Model_Type= c("Full", "NL Animal", "L0 Animal"))

## ---------------- Getting the required models and computing sample sizes ----

# A small function to select the correct model
select_model <- function(pop, type) {
    
    # Getting the relevant basename of the file
    if (type == "Full") {
        if (pop == "Moulis") {
            base <- "Moulis_comparison_LOO_model_log"
        } else if (pop == "ROB") {
            base <- "model_growth_final"
        }
    } else if (type == "NL Animal") {
        if (pop == "Moulis") {
            base <- "Moulis_animal_model_nonlinear"
        } else if (pop == "ROB") {
            base <- "animal_model_nonlinear"
        }
    } else if (type == "L0 Animal") {
        if (pop == "Moulis") {
            base <- "Moulis_animal_model_L0"
        } else if (pop == "ROB") {
            base <- "animal_model_L0_1sex"
        }
    }
    
    # Reading and returning the model
    readRDS(here("Output", str_c(base, ".rds")))
}

# Getting the corresponding model
tbl_samp <-
    tbl_samp |>
    mutate(Model = map2(Population, Model_Type,
                        select_model,
                        .progress = TRUE),
           N     = map_int(Model, \(mod_) nrow(mod_[["data"]])),
           N_ID  = map_int(Model, \(mod_) with(mod_[["data"]], n_distinct(ID))),
           Model = NULL)

## -------------------------------------------- Generating the LaTeX table ----

# Some formatting
tbl_samp <-
    tbl_samp |>
    mutate(Population = factor(Population, levels = c("ROB", "Moulis")),
           Model_Type = factor(Model_Type, levels = c("Full", "L0 Animal", "NL Animal")) |>
                        fct_recode(`Full (model selection)`  = "Full",
                                   `Animal model ($L_0$)` = "L0 Animal",
                                   `Animal model ($k$ and $L_{\\infty}$)` = "NL Animal")) |>
    arrange(Population, Model_Type) |>
    pivot_wider(names_from = Population,
                values_from = starts_with("N"),
                names_vary = "slowest")

# Generating the LaTeX table
tex_samp <-
    tbl_samp |>
    kable(format        = "latex",
          booktabs      = TRUE,
          digits        = 1,
          escape        = FALSE,
          col.names     = colnames(tbl_samp) |>
                          str_remove("Model_Type") |>
                          str_replace("N", "$N$") |>
                          str_replace("_ID", " Ind.") |>
                          str_remove("_(ROB|Moulis)")) |>
    column_spec(1, bold = TRUE, background = "blue!5!white!90!black") |>
    row_spec(0, bold = FALSE, background = "blue!15!gray") |>
    add_header_above(c("", "Wild" = 2, "Mesocosm" = 2),
                     background = "blue!15!gray",
                     bold = TRUE,
                     line = FALSE,
                     line_sep = 0) |>
    kable_styling() |>
    str_remove_all("\\\\(begin|end)\\{table\\}(\n)?")

write_file(tex_samp, file = here("Tables/Sample_sizes.tex"))
