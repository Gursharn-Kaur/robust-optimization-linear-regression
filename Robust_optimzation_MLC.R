rm(list=ls())
library(zoo)
library(dplyr)
library(data.table)
library(tidyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(tidycensus)

input_data_path <- "/Users/gursharnkaur/Documents/UVA/projects/WW_Optimization/Optimal-Wastewater-Sampling/Experiments/VA_FIPS_hospitalizations/input_data/"
output_data_path <- "/Users/gursharnkaur/Documents/UVA/projects/WW_Optimization/Optimal-Wastewater-Sampling/Experiments/VA_FIPS_hospitalizations/output_data/"


AR_residuals <- function(data, ar.order = 1){
  y <- data[, hosp]
  
  AR_model <- arima(y, order= c(ar.order, 0, 0))  # Fit an AR model
  y_resd <- AR_model$resid # AR residuals
  
  # Append AR results to the data frame
  AR_model_Rsq = 1- sum( (AR_model$resid[-1])^2)/ sum( (y-mean(y))^2) 
  AR_model_AdjRsq = 1-(1-AR_model_Rsq)*(AR_model$nobs-1)/ (AR_model$nobs- length(AR_model$coef) -1)
  
  return(list(y_resd = y_resd, AR_model_Rsq = AR_model_Rsq, AR_model_AdjRsq=AR_model_AdjRsq  ))
}


R2_subset <- function(dt, S){
  if(length(S) == 0){
    Rsq <- AR_residuals(dt)$AR_model_Rsq 
    Rsq_resd <- 0} else{
      dt[, hosp_AR_resd := AR_residuals(dt)$y_resd]
      formula <- as.formula(paste("hosp_AR_resd~", paste(S, collapse = " + ")))
      lm_model <-  lm(formula, data = dt )
      TSS <- sum( ((dt[,hosp] - mean(dt[, hosp]))^2))
      TSS_AR_resd <- sum( ((dt[,hosp_AR_resd] - mean(dt[, hosp_AR_resd]))^2))
      Rsq <- 1 - sum( (lm_model$residuals)^2)/ TSS
      Rsq_resd <-  1 - sum( (lm_model$residuals)^2)/ TSS_AR_resd
    }
  return(list(Rsq = Rsq, Rsq_resd = Rsq_resd))
}


f <- list(
  function(S) R2_subset(data[["total"]], S)$Rsq_resd,
  function(S) R2_subset(data[["adult"]], S)$Rsq_resd,
  function(S) R2_subset(data[["pediatric"]], S)$Rsq_resd
)

# Function to compute f_C
f_C <- function(S, C, f_list) {
  truncated_values <- sapply(f_list, function(func) min(func(S), C))
  return(mean(truncated_values))
}


# Implement Algorithm 1 (MLC Set Cover Algorithm)
mlc_set_cover_algorithm <- function(V, g = f_C, A, B, lambda, epsilon=.3) {
  # Initialize
  S_total <- list()  # List to store the final subsets
  S_current <- integer(0)  # Current subset
  i <- 1
  
  while (g(S_current, C, f) < g(V, C, f) - epsilon & all(A %*% as.integer(V %in% S_current) <= B) ) {
    # Define gi as the marginal gain function
    gi <- function(S) {
      return(g(union(S_current, S), C, f) - g(S_current, C, f))
    }
    
    
    # Initialize Si and weights
    Si <- integer(0)
    w <- 1 / B
    
    while (sum(B * w) <= lambda && length(S_current) < length(V)) {
      # Find the best element to add
      
      j <- which.min(sapply(V[!(V %in% Si)], function(j) {
        sum(A[, j] * w) / (gi(c(S_current, j)) - gi(S_current))
      }))
      
      # Update Si and weights
      Si <- union(Si, names(j))
      w <- w * lambda^( A[,j ] / B)
    }
    
    Si_binary <- as.integer(V %in% Si)
    
    if (all(A %*% Si_binary <= B) | length(Si) == 1 ) {
      Si <- Si
    } else {
      # Constraint check and adjustment
      if (gi(setdiff(Si, predictors[j])) >= gi(predictors[j])) {
        Si <- setdiff(Si, predictors[j])
      } else {
        Si <- predictors[j]
      }
    }
    
    # Update the total solution
    S_total <- union(S_total, list(Si))
    S_current <- union(S_current, Si)
    i <- i + 1
  }
  return(S_total)
}



robust_MLC_set_cover <- function(V, delta = .2){
  
  S_best = vector(); C_min = 0; C_max = 1
  while( (C_max - C_min) > delta){
    C = (C_max + C_min)/2
    
    MLC_cover_output <-  mlc_set_cover_algorithm(V, g = f_C, A = A, B = B, lambda = lambda, epsilon = epsilon)
    
    S_binary <- as.integer(V %in% unlist(MLC_cover_output))
    
    if (all(A %*% S_binary > length(MLC_cover_output)* B) ) {
      C_max <- C }else {
        C_min <- C; S_best <- unlist(MLC_cover_output)
      }
    #print(S_best)
  }
  return(S_best)
}



##### Implementation

#Input time series data 
scenario_set <- c("adult", "pediatric", "total")
data <- list()
for( sc in scenario_set){
  data[[sc]] <- fread(paste0(input_data_path, "data_", sc, ".csv"))
}



# Get predictor names starting with X_
predictors <- grep("^X_", colnames(data[["total"]]), value = TRUE)

# Define constraints 
n <- length(predictors); m <- 2

A = matrix(0, ncol = n, nrow = m) # constraint matrix 
colnames(A) <- predictors

A[1,] <- 1 # first constraint 

#population data required for second constraint

My_key <- "396e95bbfe59ba3816b33d8dac2390715b5ee429"

va_fips_population <- get_acs(
  geography = "county",
  variables = "B01003_001E",  # Total population
  state = "51",               # FIPS code for Virginia
  year = 2019,                # Year of data
  survey = "acs5",             # American Community Survey 5-Year Estimates,
  key = My_key
)

# Process the data
pop_fips <- va_fips_population %>%
  select(GEOID, estimate) %>%
  rename(geo_value = GEOID, population = estimate) %>%
  mutate(geo_value = paste0("X_", geo_value),
         proportion_pop = population / sum(population)) %>%
  data.table()

#define second constraint
for (predictor in predictors) {
  if (predictor %in% pop_fips[,geo_value]) {
    A[2, predictor] <- 1 - pop_fips[geo_value== predictor , proportion_pop]
  }
}

# Example inputs
V <- predictors  
A <- A  # Constraint matrix A
lambda <- 2.1 # Update factor
epsilon <- 0.1  # Threshold for convergence
scenario_set <- scenario_set


# Run the MLC set cover algorithm
result_robust_MLC = list(); Rsq_MLC <- numeric(10)
for( k in 1: 10){
  # Fix budget 
  B <- c(k, k+.1 - 1) 
  result_robust_MLC[[k]]  <- robust_MLC_set_cover(V, delta = .1)
  print(k)
}

output_table <- data.table()
for( k in 1:10){
  Rsq_values <- numeric(length(scenario_set))
  
  # Compute Rsq values for all scenarios
  for (j in 1:length(scenario_set)) {
    Rsq_values[j] <- sR2_subset(data[[scenario_set[j]]], result_robust_MLC[[k]])$Rsq
  }
  
  # Save the minimum Rsq value across all scenarios for the current k
  Rsq_MLC[k] <- min(Rsq_values)
  
  output_table <- rbind(output_table, 
                        data.table(k=k, predictors = paste0(result_robust_MLC[[k]], collapse = ", "), min_Rsq = Rsq_MLC[k])
  )
}


kable(output_table , digits=4)





