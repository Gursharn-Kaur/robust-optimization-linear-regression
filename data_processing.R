rm(list=ls())
library(zoo)
library(dplyr)
library(data.table)
library(tidyr)
library(reshape2)


#Input time series data 
scenario_set <- c("adult", "pediatric", "total")

data_hosp <- fread("/Users/gursharnkaur/Downloads/COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_VAcounty.csv")

setnames(data_hosp, old = c("date",  "fips_code",                                            
                            "previous_day_admission_adult_covid_confirmed_7_day_sum",    
                            "previous_day_admission_pediatric_covid_confirmed_7_day_sum",
                            "previous_day_admission_total_covid_confirmed_7_day_sum") ,
         new = c("time_value", "geo_value", "adult", "pediatric", "total"),
         skip_absent = TRUE)

data <- list()
data[["adult"]] <- select(data_hosp, c("time_value", "geo_value", "adult"))
data[["pediatric"]] <- select(data_hosp, c("time_value", "geo_value", "pediatric"))
data[["total"]] <- select(data_hosp, c("time_value", "geo_value", "total"))

for( sc in scenario_set){
  filtered_geo <- data[[sc]] %>%
    count(geo_value) %>%
    filter(n >= 180) %>%
    pull(geo_value)
  
  filtered_data <- data[[sc]] %>%
    filter(geo_value %in% filtered_geo)
  
  data[[sc]] <- filtered_data %>%
    pivot_wider(
      names_from = geo_value,
      values_from = sc,
      names_prefix = "X_"
    )%>%
    data.table()
  data[[sc]] <- data[[sc]] %>%
    mutate(across(-time_value, ~ na.locf(.x, na.rm = FALSE)))
  
  data[[sc]] <- na.omit(data[[sc]])
  
  data[[sc]] <- data[[sc]] %>%
    mutate(hosp = rowSums(select(., starts_with("X_")), na.rm = TRUE))
}

write.csv(data[["adult"]], 
          "/Users/gursharnkaur/Documents/UVA/projects/WW_Optimization/Optimal-Wastewater-Sampling/Experiments/VA_FIPS_hospitalizations/input_data/data_adult.csv",
          row.names = FALSE)

write.csv(data[["pediatric"]], 
          "/Users/gursharnkaur/Documents/UVA/projects/WW_Optimization/Optimal-Wastewater-Sampling/Experiments/VA_FIPS_hospitalizations/input_data/data_pediatric.csv",
          row.names = FALSE)
write.csv(data[["total"]], 
          "/Users/gursharnkaur/Documents/UVA/projects/WW_Optimization/Optimal-Wastewater-Sampling/Experiments/VA_FIPS_hospitalizations/input_data/data_total.csv",
          row.names = FALSE)
