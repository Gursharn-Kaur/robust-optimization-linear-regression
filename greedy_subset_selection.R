# -------------------------------------------------------------------
#  Robust subset-selection demo (anonymous version)
# -------------------------------------------------------------------

rm(list = ls())

# --------------------------  Libraries  ----------------------------
library(zoo)
library(dplyr)
library(data.table)
library(tidyr)
library(leaps)
library(knitr)
library(reshape2)
library(ggplot2)
library(plotly)
# (covidcast removed; not used below.  Add back if required.)

# ---------------------  User-configurable paths  -------------------
#  !!  EDIT ONLY THESE TWO LINES  !!
input_data_path  <- file.path("input_data")       # folder with data_*.csv
output_data_path <- file.path("output_results")   # folder for outputs
dir.create(output_data_path, showWarnings = FALSE, recursive = TRUE)

# -----------------------  Helper functions  ------------------------
calculate_R2 <- function(predictor_subset, y_resd, data, TSS) {
  f  <- as.formula(paste("y_resd ~", paste(predictor_subset, collapse = " + ")))
  lm_model <- lm(f, data = data)
  n  <- length(lm_model$residuals)
  Rsq <- 1 - sum(lm_model$residuals^2) / TSS
  Adj_Rsq <- 1 - (1 - Rsq) * (n - 1) / (n - length(predictor_subset) - 1)
  Rsq_resd <- summary(lm_model)$r.squared
  
  coef_tbl <- broom::tidy(lm_model)
  signif_predictors <- coef_tbl$term[coef_tbl$p.value < 0.05 &
                                       coef_tbl$term %in% predictor_subset]
  
  data.table(
    k  = length(predictor_subset),
    predictors             = paste(predictor_subset, collapse = ", "),
    significant_predictors = paste(signif_predictors, collapse = ", "),
    Rsq        = Rsq,
    Rsq_resd   = Rsq_resd,
    Adj_Rsq    = Adj_Rsq
  )
}

AR_residuals <- function(data, ar.order = 1) {
  y <- data[ , hosp]
  ar_fit <- arima(y, order = c(ar.order, 0, 0))
  resid  <- ar_fit$resid
  Rsq    <- 1 - sum(resid[-1]^2) / sum((y - mean(y))^2)
  Adj_Rsq <- 1 - (1 - Rsq) * (ar_fit$nobs - 1) /
    (ar_fit$nobs - length(ar_fit$coef) - 1)
  list(y_resd = resid, AR_model_Rsq = Rsq, AR_model_AdjRsq = Adj_Rsq)
}

# --------------------  Greedy-subset algorithm  --------------------
greedy_algorithm <- function(dt, lag = 1) {
  setorder(dt, time_value)
  y <- dt[ , hosp]
  
  dt_lag <- dt
  for (l in seq_len(lag))
    dt_lag <- dt_lag %>% mutate(across(-time_value, dplyr::lag))
  
  dt_all_var <- dt_lag[ , !c("time_value", "hosp"), with = FALSE]
  predictors <- names(dt_all_var)
  
  ar_info <- AR_residuals(dt_lag)
  y_resd  <- ar_info$y_resd
  TSS     <- sum((y - mean(y))^2)
  
  result_tbl <- data.table(
    k        = 0,
    predictors             = "hosp",
    significant_predictors = "hosp",
    Rsq      = ar_info$AR_model_Rsq,
    Rsq_resd = 0
  )
  
  best_subset <- character(0)
  best_Rsq    <- 0
  
  for (k in seq_along(predictors)) {
    if (k == 1) {
      subsets <- as.list(predictors)
    } else {
      remain <- setdiff(predictors, best_subset)
      subsets <- lapply(remain, function(x) c(best_subset, x))
    }
    
    for (sub in subsets) {
      obj <- calculate_R2(sub, y_resd, dt_all_var, TSS)
      if (obj$Rsq > best_Rsq) {
        best_Rsq    <- obj$Rsq
        best_subset <- sub
        best_signif <- obj$significant_predictors
        best_Rsq_resd <- obj$Rsq_resd
      }
    }
    
    result_tbl <- rbind(
      result_tbl,
      data.table(
        k        = length(best_subset),
        predictors             = paste(best_subset, collapse = ", "),
        significant_predictors = best_signif,
        Rsq      = best_Rsq,
        Rsq_resd = best_Rsq_resd
      )
    )
  }
  result_tbl
}

# ---------------  Load scenario data (adult/pediatric/total) -------
scenario_set <- c("adult", "pediatric", "total")
data_list <- lapply(
  scenario_set,
  function(sc) fread(file.path(input_data_path, paste0("data_", sc, ".csv")))
)
names(data_list) <- scenario_set

# -------------------------  Run greedy search  ---------------------
for (sc in scenario_set) {
  res <- greedy_algorithm(data_list[[sc]])
  fwrite(
    res,
    file = file.path(output_data_path,
                     paste0("greedy_results_", sc, ".csv"))
  )
}

# ----------------------------  Done  -------------------------------
message("All scenario runs completed.  Results saved in: ", output_data_path)
