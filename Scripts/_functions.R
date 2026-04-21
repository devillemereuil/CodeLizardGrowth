##########################################################################################
##                                 Functions for the                                    ##
##                             "Genetics of growth analysis"                            ##
##                          Anaïs Aragon & Pierre de Villemereuil                       ##
##                                       (2023)                                         ##
##########################################################################################

## ------------------------------------------------------- Maths functions ----

## Logistic growth function
# Equation of the model (Schoener et al. 1978)
# 3 parameters :
#   L_0 : size at birth
#   L_max : maximal size
#   k > 0 : growth rate
f_log <- function(x, L0, Lmax, k) {
    Lmax  / (1 + ((Lmax - L0) / L0) * exp(- k * x))
}

## von Bertalanffy growth function
# Equation of the model (Schoener et al. 1978)
# 3 parameters :
#   L_0 : size at birth
#   L_max : maximal size
#   k > 0 : growth rate
f_vb <- function(x, L0, Lmax, k) {
    Lmax - (Lmax - L0) * exp(- k * x)
}

## Gompertz growth function
# Equation of the model (Chen Yang et al. 2018)
# 3 parameters :
#   L_0 : size at birth
#   L_max : maximal size
#   k > 0 : growth rate
f_gz <- function(x, L0, Lmax, k) {
    Lmax * exp(log(L0 / Lmax) * exp(- k * x))
}

## ------------------------------------------------------- Utils functions ----

## Function to collapse a vector with a unique level to a unique value
# Args : - vec: a vector with possible a unique level (e.g. c(1, 1, 1))
# Value : a unique value to replace the vector (if more than 1, return NA)
collapse_to_one <- function(vec) {
    n <- n_distinct(vec, na.rm = TRUE)
    n_na <- sum(is.na(vec))
    if (n == 1 & n_na == 0) {
        return(unique(vec))
    } else if (n == 0) {
        return(NA)
    } else {
        return(
            sort(table(vec),decreasing=TRUE)[1] |>
                names()
        )
    }
}

## Convert to character with significant digits and trailing zeros
# Args : - num: a number
#        - digits: significant figures
# Value : a character string with significant digits (and trailing zeros)
to_signif <- function(num, digits) {
    sprintf(paste0('%#.', digits, 'g'), signif(num, digits))
}

## LaTeX scientific notation
# Args : - num: a number
#        - digits: significant figures
# Value : a character string with LaTeX formatted scientific notations
latex_sci <- function(num, digits) {
    exp     <- floor(log10(num))
    base    <- to_signif(num * 10^(-exp), digits)
    stringr::str_c(base,"\\text{\\textsc{e",exp,"}}")
}


## -------------------------------------------- Models wrangling functions ----

## Extract sex-specific intercepts from a model
# Args : - mod: a brmsfit model object from this project
#        - params: a list of the non-linear parameters (default NULL for linear models)
# Value : the MCMC output of the sex-specific intercepts
get_sex_intercept <- function(mod, params = NULL, corr = FALSE) {

    # Extracting the fixed effects
    fe <- fixef(mod, summary = FALSE)

    # Computing the sex-specific intercepts
    if (is.null(params)) {

        # Extract on intercept- and sex-related effects
        fe <- fe[ , c("Intercept", "SexM")]

        # Compute the male intercept
        fe[ , "SexM"] <- fe[ , "Intercept"] + fe[ , "SexM"]
        colnames(fe) <- c("Int_SexF", "Int_SexM")

        out <- fe

    } else { # If the model is non-linear, working parameter by parameter

        # Check whether Sex is the sole fixed effect
        if (ncol(fe) != 2 * length(params)) {
            stop("More than two fixed-effect parameters by parameter,
                 is there another fixed effect than Sex?")
        }

        # Now working parameter-by-parameter to compute male intercept
        for (p in params) {
            fe[ , paste0(p, "_SexM")] <-
                fe[ , paste0(p, "_Intercept")] + fe[ , paste0(p, "_SexM")]
        }
        colnames(fe) <-
            stringr::str_replace(colnames(fe),
                                 "Intercept",
                                 "Int_SexF")
        colnames(fe) <-
            stringr::str_replace(colnames(fe),
                                 "SexM",
                                 "Int_SexM")

        list_F <- paste0(params, "_Int_SexF")
        list_M <- paste0(params, "_Int_SexM")

        fe <- as_tibble(fe)

        if (corr) {
            # Packing parameters together by sex
            out <-
                tibble(Int_SexF = pmap(fe[list_F], c) |>
                            map(\(vec) { names(vec) <- params; return(vec) }),
                    Int_SexM = pmap(fe[list_M], c) |>
                            map(\(vec) { names(vec) <- params; return(vec) }))
        } else {
            out <-
                mutate(fe, Iter = 1:n()) |>
                pivot_longer(-Iter,
                             names_to = c("Parameter", ".value"),
                             names_pattern = "^([^_]*)_(.*)$")
        }
    }

    return(as_tibble(out))
}

## Extract the variance components from a model
# Args : - mod: a brmsfit model object from this project
#        - params: a list of the non-linear parameters (default NULL for linear models)
# Value : the MCMC output of the random effects variance components
get_var_random <- function(mod, params = NULL, corr = FALSE) {

    ## Extracting the fixed effects
    vars <- VarCorr(mod, summary = FALSE)

    # Computing the sex-specific intercepts
    if (is.null(params)) {

        ## Get random effect variances
        out <-
            vars |>
            purrr::map("sd") |>
            purrr::imap(\(df, name) {
                tibble::as_tibble(df) |>
                    dplyr::rename_with(\(col_) stringr::str_replace(col_, "V1", "Intercept")) |>
                    dplyr::mutate(Iter = 1:n()) |>
                    tidyr::pivot_longer(-Iter,
                                        names_to = "Sex",
                                        values_to = "Estimate") |>
                    dplyr::mutate(Estimate = Estimate^2,
                                  Variance = dplyr::case_match(name,
                                                               "ID"          ~ "V_A",
                                                               "ID_Mom"      ~ "V_M",
                                                               "Year_fac"    ~ "V_Y",
                                                               "residual__"  ~ "V_R"))
            }) |>
            purrr::list_rbind() |>
            tidyr::pivot_wider(names_from = c(Variance, Sex),
                               values_from = "Estimate") |>
            dplyr::rename_with(\(col_) stringr::str_remove(col_, "_Intercept"))

        ## Get fixed effect variances
        # Fixed effects (without intercept and sex effect)
        fe <- fixef(mod, summary = FALSE) |>
              as.data.frame() |>
              select(-Intercept, -SexM, -starts_with("sigma"))
        # Extracting a design matrix from the data
        # NOTE: Will not work for dummy variables but brms does not offer a simple
        #       way to compute the design matrix
        X <- mod[["data"]] |>
             dplyr::select(colnames(fe)) |>
             dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
             as.matrix()
        # To account for the uncertainty, we need the SE^2 matrix
        S <- vcov(mod)[colnames(fe), colnames(fe), drop = FALSE]
        var_uncert <- sum(cov(X) %*% S)
        # Now we can compute the fixed effects
        vf <- apply(fe, 1, \(beta) { var(X %*% beta) }) - var_uncert

        out <- dplyr::mutate(out, V_F = vf, .before = 1)

        ## Adding residual variance if missing
        if (!("V_R" %in% colnames(out))) {
            sigma <-
                fixef(mod, summary = FALSE) |>
                as.data.frame() |>
                select(starts_with("sigma")) |>
                transmute(V_R_SexF = exp(sigma_SexF)^2,
                          V_R_SexM = exp(sigma_SexM)^2)

            out <-
                bind_cols(out, sigma)
        }

    } else if (corr) { # If the model is non-linear, working parameter by parameter
        # Removing the residual variance
        vars[["residual__"]] <- NULL

        # Generating the cross-values for parameters
        cross <-
            tidyr::crossing(P1 = params, P2 = params) |>
            dplyr::filter(P1 != P2) |>
            with(paste(P1, P2, sep = "."))
        full <-
            tidyr::crossing(P1 = c(params), P2 = c(params)) |>
            with(paste(P1, P2, sep = "."))

        ## Get random effect variances
        out <-
            vars |>
            purrr::imap(\(lst, name) {
                if (is.null(lst[["cov"]])) {
                    tbl <-
                        lst[["sd"]] |>
                        tibble::as_tibble() |>
                        dplyr::mutate(across(everything(), \(x) { x^2 })) |>
                        dplyr::rename_with(\(chr) { paste0(chr,".",chr) })
                    tbl[ , cross] <- 0

                } else {
                    tbl <-
                        lst[["cov"]] |>
                        tibble::as_tibble()
                }

                tbl <-
                    tbl |>
                    dplyr::rename_with(\(chr) gsub("_Intercept", "", chr))

                out <- tibble({{name}} := map(1:nrow(tbl),
                                   \(i) { out <- matrix(unlist(tbl[i, full]),
                                                 ncol = length(params),
                                                 nrow = length(params));
                                          colnames(out) <-
                                              rownames(out) <-
                                                  params;
                                          return(out)}))
                return(out)
            }) |>
            purrr::list_cbind() |>
            tidyr::unpack(everything()) |>
            dplyr::rename_with(\(vec) {
                dplyr::case_match(vec,
                                  "ID"          ~ "V_A",
                                  "ID_Mom"      ~ "V_M",
                                  "ID_PE"       ~ "V_PE",
                                  "Year_fac"    ~ "V_Y")

            })

        # NOTE Vf is not computed here, it's not needed in our particular case and
        # would be a lot of work
    } else {
        # Removing the residual variance
        vars[["residual__"]] <- NULL

        out <-
            vars |>
            purrr::map("sd") |>
            purrr::imap(\(df, name) {
                tibble::as_tibble(df) |>
                    dplyr::mutate(Iter = 1:n()) |>
                    tidyr::pivot_longer(-Iter,
                                        names_to = "Parameter",
                                        values_to = "Estimate") |>
                    separate(col = "Parameter", sep = "_", into = c("Parameter", "Level")) |>
                    dplyr::mutate(Estimate = Estimate^2,
                                  Variance = dplyr::case_match(name,
                                                               "ID"          ~ "V_A",
                                                               "ID_Mom"      ~ "V_M",
                                                               "ID_PE"       ~ "V_PE",
                                                               "Year_fac"    ~ "V_Y"))
            }) |>
            purrr::list_rbind() |>
            tidyr::pivot_wider(names_from = c(Variance, Level),
                               values_from = "Estimate") |>
            dplyr::rename_with(\(col_) stringr::str_remove(col_, "_Intercept"))
    }

    return(as_tibble(out))
}

## Add extra parameters from the intercept and variance output
# Args : - out: a tibble containing the combined output of get_sex_intercept()
#               and get_var_random()
#        - params: a list of the non-linear parameters (default NULL for linear models)
# Value : the MCMC output of the random effects variance components
add_deriv_parameters <- function(out, params = NULL, corr = FALSE) {

    # Computing the sex-specific intercepts
    if (is.null(params) | !corr) {
        if ("Sex" %in% colnames(out)) {
            out[["Mean"]] <- out[["Int"]]
        } else {
            out[["Mean"]] <-
                dplyr::select(out, Int_SexF, Int_SexM) |>
                rowMeans()
        }
        out[["V_P"]] <-
            dplyr::select(out, starts_with("V_")) |>
            rowSums()
        out[["Herit"]] <- with(out, V_A / V_P)
        out[["Evolv"]] <- with(out, V_A / (Mean^2))
    } else { # If the model is non-linear, working parameter by parameter
        if ("Sex" %in% colnames(out)) {
            out[["Mean"]] <- out[["Int"]]
        } else {
            out[["Mean"]] <-
                purrr::map2(out[["Int_SexF"]], out[["Int_SexM"]],
                            \(vec_f, vec_m) { 0.5 * (vec_f + vec_m) })
        }
        out[["V_P"]] <-
            purrr::pmap(out |> select(starts_with("V")),
                        \(...) { purrr::reduce(list(...), `+`) })
        out[["Herit"]] <-
            purrr::map2(out[["V_A"]], out[["V_P"]],
                        \(G, P) { h2 <- diag(G) / diag(P); names(h2) <- params; h2 })
        out[["Evolv"]] <-
            purrr::map2(out[["V_A"]], out[["Mean"]],
                        \(G, mu) { I <- diag(G) / (mu^2); names(I) <- params; I })
    }

    return(as_tibble(out))
}

## Summarise HMC chains
# Args : - vec: a vector containing the iterations for a specific parameter
#        - with_p: should a p-value against H0 = 0 be computed?
# Value : a tbl summarising the chains
summarise_chains <- function(vec, with_p = FALSE) {
    out <-
        tibble::tibble(Mean         = mean(vec),
                       Median       = median(vec),
                       Low          = quantile(vec, p = 0.025),
                       Up           = quantile(vec, p = 0.975))

    if (with_p) {
        out <-
            out |>
            dplyr::mutate(P_val     = 2 * sum(sign(vec) != sign(Mean)) / length(vec),
                          P_val_min = 1/length(vec))
    }

    return(out)
}

## ----------------------------------------------------- Output formatting ----

## Group columns for add_header_above() in kableExtra
# Args : - names: a vector of colnames
# Value : a named vector counting the subsequent identical columns to pass to add_header_above()
group_col <- function(names) {
    names <- str_replace(names, "^$", " ")
    grp <- rle(names)
    return(setNames(grp[["lengths"]], grp[["values"]]))
}
