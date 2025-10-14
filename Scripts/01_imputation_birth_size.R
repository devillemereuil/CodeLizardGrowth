conflicted::conflicts_prefer(dplyr::filter)

###############################################################################
##                Imputation of the size at birth when missing               ##
##                       "Genetics of growth analysis"                       ##
##                     Anaïs Aragon & Pierre de Villemereuil                 ##
##                                  (2023)                                   ##
###############################################################################

# Packages required
library(tidyverse)
library(here)
library(mice)

# Setting up a (lucky) seed
seed <- 777
set.seed(seed)

# Source the custom functions for this analysis
source("_functions.R")

## --------------------------------------- Loading and formatting the data ----

# Loading the datasets
tbl_capt <- readRDS(here("Data/capture_clean.rds"))
tbl_repro <- readRDS(here("Data/repro_clean.rds"))

# Setting up the levels for the years 
lev_years <- (1993:2019) |> as.character()


## ---------------------------------------- Imputing missing data with mice ---

# Generating a wide table with NB, JC, SA size for each ID
tbl_size <-
    tbl_capt |>
    arrange(Date_Capt, ID) |>
    summarise(ID         = unique(ID),
              Sex        = unique(Sex),
              Year_Birth = unique(Year - Age),
              Size_NB    = first(Size[which(Stage == "NB")]),
              Size_JC    = first(Size[which(Stage == "JC")]),
              Month_SA   = Date_Capt[which(Stage == "SA")] |>
                                       month() |>
                                       first(),
              Size_SA    = first(Size[which(Stage == "SA")]),
              .by = ID)

# Imputing the data using chained equations and mice
imput <- mice(tbl_size, m = 20, seed = seed)

# Extracting the lines where Size_NB was imputed
rows_imput <- imput[["where"]][ , "Size_NB"]

# Getting the average imputations
tbl_avg_imput <-
    imput |>
    complete(action = "long") |>
    as_tibble() |>
    summarise(ID        = unique(ID),
              Sex       = unique(Sex),
              Stage     = "NB",
              Age       = 0,
              Year      = unique(Year_Birth),
              Size      = mean(Size_NB),
              Size_SE   = sd(Size_NB),
              .by = .id) |>
    select(-.id)

# Getting a list of imputed Size_NB tbls
lst_imput <-
    imput |>
    complete(action = "all") |>
    map(as_tibble) |>
    map(\(df_) {
        df_ |>
        transmute(ID    = ID,
                  Sex   = Sex,
                  Stage = "NB",
                  Age   = 0,
                  Year  = Year_Birth,
                  Size  = Size_NB) |>
        filter(rows_imput) |>
        mutate(Imputed_Size_NB = TRUE)
    })

## -------------------------------------- Merging with the Capture dataset ----

# Dataset with average predictions
tbl_capt_avgimp <-
    bind_rows(
        tbl_capt |> mutate(Size_SE = 0),
        tbl_avg_imput |>
            filter(Size_SE > 0)
    ) |>
    left_join(tbl_avg_imput |> select(ID, Size_NB = Size, Size_NB_SE = Size_SE)) |>
    arrange(ID, Year, Date_Capt)

# List of multiple imputations
lst_capt_imp <-
    lst_imput |>
    map(\(df_) {
        bind_rows(
            tbl_capt |> mutate(Imputed_Size_NB = FALSE),
            df_
        ) |>
        mutate(Size_NB = Size[Stage == "NB"],
               Imputed_Size_NB = any(Imputed_Size_NB),
               .by = ID) |>
        arrange(ID, Year, Date_Capt)
    })

## ------------------------------- Exporting the dataset with imputed data ----

saveRDS(tbl_capt_avgimp, file = here("Data/capture_avg_imputed.rds"))
saveRDS(lst_capt_imp, file = here("Data/capture_list_imputed.rds"))

