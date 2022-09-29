### IO problem set - warmup questions


library(tidyverse)
library(fixest)
library(modelsummary)

# load data

OTC_instruments <- read.csv('OTCDataInstruments.csv', sep = "")
OTC_headache <-read.delim("~/GitHub/PhD_IO/pset1/OTC_Headache_upd.csv", header=TRUE)
store_dem <- read.csv('store_demographics.csv')


# Reshaping so that each row corresponds to a week-store-brand-size combination

product_data <- OTC_headache %>%
  pivot_longer(sales_1:prom_11, 
               names_to = c(".value", "product-brand"),
               names_pattern =  "([A-Za-z]+)[^0-9]+([0-9]+)")

product_data <- product_data %>%
  left_join(store_dem, by = c('store'='STORE'))
