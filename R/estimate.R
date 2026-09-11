# Shared argument checking for the P_Z / Q_Z / V_Z estimators.
check_levels <- function(data, instrument, treatment,
                         instrument_levels, treatments) {
  for (col in c(instrument, treatment)) {
    if (!col %in% names(data)) {
      stop("Column ", sQuote(col), " is not present in `data`.", call. = FALSE)
    }
  }
  if (is.null(instrument_levels)) {
    instrument_levels <- sort(unique(data[[instrument]]))
  }
  if (is.null(treatments)) {
    treatments <- sort(unique(as.character(data[[treatment]])))
  }
  missing_z <- setdiff(as.character(instrument_levels),
                       as.character(unique(data[[instrument]])))
  if (length(missing_z)) {
    stop("No observations for instrument value(s): ",
         paste(missing_z, collapse = ", "),
         ". Every instrument value must be observed, because the rows of the ",
         "returned matrix have to line up with the rows of the response ",
         "matrix.", call. = FALSE)
  }
  list(instrument_levels = as.character(instrument_levels),
       treatments = as.character(treatments))
}

# Turn long marginaleffects output into an instrument-by-treatment matrix.
widen_effects <- function(long, z_col, t_col, value_col,
                          instrument_levels, treatments) {
  long <- as.data.frame(long)
  out <- matrix(NA_real_, nrow = length(instrument_levels),
                ncol = length(treatments),
                dimnames = list(instrument_levels, treatments))
  i <- match(as.character(long[[z_col]]), instrument_levels)
  j <- match(as.character(long[[t_col]]), treatments)
  keep <- !is.na(i) & !is.na(j)
  out[cbind(i[keep], j[keep])] <- long[[value_col]][keep]
  out
}

#' Estimate the conditional treatment probabilities P_Z
#'
#' `P_Z(t)[i] = P(T = t | Z = z_i)`, the probability of receiving each
#' treatment under each value of the instrument.
#'
#' @param data A data frame with one row per observation.
#' @param instrument Name of the instrument column, as a string.
#' @param treatment Name of the treatment column, as a string.
#' @param covariates Character vector of columns to adjust for. Only available
#'   when `parametric = TRUE`; supplying them otherwise is an error rather than
#'   being silently ignored.
#' @param parametric Estimate by multinomial regression rather than by cell
#'   proportions. Required for covariate adjustment.
#' @param instrument_levels,treatments Optional character vectors fixing the
#'   row and column order. Default to the sorted unique values found in `data`.
#'   Pass them explicitly to guarantee the rows line up with the rows of
#'   [response_matrix()].
#'
#' @return A numeric matrix with one row per instrument value and one column
#'   per treatment. Rows sum to one.
#'
#' @details
#' The parametric estimator fits a multinomial regression of treatment on the
#' instrument and the covariates, then averages the predicted probabilities
#' over the observed covariate distribution separately at each value of the
#' instrument. That standardisation is what makes the adjusted estimator a
#' valid input to the identification results.
#'
#' Called `MakeP_Z()` in the paper code.
#'
#' @seealso [estimate_q_z()], [estimate_p_sigma()]
#'
#' @examples
#' d <- data.frame(
#'   z = rep(1:3, each = 40),
#'   trt = rep(c("a", "b", "a", "b", "b", "a"), each = 20)
#' )
#' estimate_p_z(d, "z", "trt")
#'
#' @export
estimate_p_z <- function(data, instrument, treatment, covariates = NULL,
                         parametric = FALSE, instrument_levels = NULL,
                         treatments = NULL) {
  if (!is.null(covariates) && !parametric) {
    stop("Covariate adjustment requires `parametric = TRUE`.", call. = FALSE)
  }
  lv <- check_levels(data, instrument, treatment, instrument_levels, treatments)

  if (!parametric) {
    counts <- table(factor(as.character(data[[instrument]]),
                           levels = lv$instrument_levels),
                    factor(as.character(data[[treatment]]),
                           levels = lv$treatments))
    out <- unclass(prop.table(counts, margin = 1L))
    out[is.na(out)] <- 0
    dimnames(out) <- list(lv$instrument_levels, lv$treatments)
    return(out)
  }

  fit_data <- data
  fit_data[[instrument]] <- factor(as.character(fit_data[[instrument]]),
                                   levels = lv$instrument_levels)
  fit_data[[treatment]] <- factor(as.character(fit_data[[treatment]]),
                                  levels = lv$treatments)
  form <- stats::as.formula(
    paste(treatment, "~", paste(c(instrument, covariates), collapse = " + "))
  )
  model <- nnet::multinom(form, data = fit_data, trace = FALSE)
  long <- marginaleffects::avg_predictions(model, variables = instrument)
  out <- widen_effects(long, instrument, "group", "estimate",
                       lv$instrument_levels, lv$treatments)
  out[is.na(out)] <- 0
  out
}

#' Estimate the outcome-probability products Q_Z
#'
#' `Q_Z(t)[i] = E(Y | Z = z_i, T = t) * P(T = t | Z = z_i)`. Ratios of `Q_Z`
#' to `P_Z`, both projected through `b B_t^+`, give the local average
#' treatment responses.
#'
#' @inheritParams estimate_p_z
#' @param outcome Name of the outcome column, as a string.
#' @param family Family passed to [stats::glm()] for the parametric outcome
#'   model, for example `"binomial"` for a binary outcome.
#' @param p_z A matrix from [estimate_p_z()]. Required when
#'   `parametric = TRUE`, because the adjusted outcome means are multiplied by
#'   the adjusted treatment probabilities.
#'
#' @return A numeric matrix with one row per instrument value and one column
#'   per treatment.
#'
#' @details
#' Without adjustment this reduces to the sum of the outcome within each
#' instrument-by-treatment cell, divided by the number of observations at that
#' value of the instrument. Cells with no observations contribute zero, which
#' is right here because `P(T = t | Z = z)` is zero for them too.
#'
#' Called `MakeQ_Z()` in the paper code.
#'
#' @export
estimate_q_z <- function(data, instrument, treatment, outcome,
                         covariates = NULL, parametric = FALSE, family = NULL,
                         p_z = NULL, instrument_levels = NULL,
                         treatments = NULL) {
  if (!is.null(covariates) && !parametric) {
    stop("Covariate adjustment requires `parametric = TRUE`.", call. = FALSE)
  }
  if (!outcome %in% names(data)) {
    stop("Column ", sQuote(outcome), " is not present in `data`.", call. = FALSE)
  }
  lv <- check_levels(data, instrument, treatment, instrument_levels, treatments)
  z_fac <- factor(as.character(data[[instrument]]),
                  levels = lv$instrument_levels)
  t_fac <- factor(as.character(data[[treatment]]), levels = lv$treatments)

  if (!parametric) {
    totals <- tapply(data[[outcome]], list(z_fac, t_fac), sum)
    totals[is.na(totals)] <- 0
    out <- totals / as.numeric(table(z_fac))
    dimnames(out) <- list(lv$instrument_levels, lv$treatments)
    return(out)
  }

  if (is.null(p_z)) {
    stop("`p_z` must be supplied when `parametric = TRUE`.", call. = FALSE)
  }
  fit_data <- data
  fit_data[[instrument]] <- z_fac
  fit_data[[treatment]] <- t_fac
  rhs <- paste0(instrument, " * ", treatment)
  if (!is.null(covariates)) rhs <- paste(c(rhs, covariates), collapse = " + ")
  form <- stats::as.formula(paste(outcome, "~", rhs))
  model <- stats::glm(form, data = fit_data, family = family)
  long <- marginaleffects::avg_predictions(
    model, variables = c(instrument, treatment)
  )
  means <- widen_effects(long, instrument, treatment, "estimate",
                         lv$instrument_levels, lv$treatments)
  means[is.na(means)] <- 0
  means * p_z[lv$instrument_levels, lv$treatments, drop = FALSE]
}

#' Estimate the conditional outcome variances V_Z
#'
#' `V_Z(t)[i] = Var(Y | Z = z_i, T = t) * P(T = t | Z = z_i)^2`, the
#' contribution of each instrument-by-treatment cell to the variance of
#' [estimate_q_z()].
#'
#' @inheritParams estimate_q_z
#'
#' @return A numeric matrix with one row per instrument value and one column
#'   per treatment. Cells holding fewer than two observations contribute zero,
#'   since no within-cell variance can be estimated from them.
#'
#' @details Called `MakeV_Z()` in the paper code. There is no parametric
#'   version.
#'
#' @export
estimate_v_z <- function(data, instrument, treatment, outcome,
                         instrument_levels = NULL, treatments = NULL) {
  if (!outcome %in% names(data)) {
    stop("Column ", sQuote(outcome), " is not present in `data`.", call. = FALSE)
  }
  lv <- check_levels(data, instrument, treatment, instrument_levels, treatments)
  z_fac <- factor(as.character(data[[instrument]]),
                  levels = lv$instrument_levels)
  t_fac <- factor(as.character(data[[treatment]]), levels = lv$treatments)

  variances <- tapply(data[[outcome]], list(z_fac, t_fac), stats::var)
  counts <- unclass(table(z_fac, t_fac))
  out <- variances * (counts / rowSums(counts))^2
  out[is.na(out)] <- 0
  dimnames(out) <- list(lv$instrument_levels, lv$treatments)
  out
}
