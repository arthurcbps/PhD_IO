### IO problem set - warmup questions


library(tidyverse)
library(fixest)
library(modelsummary)
library(janitor)
library(kableExtra)

# load data

OTC_instruments <- read.csv("OTCDataInstruments.csv", sep = "")
OTC_headache <- read.delim("~/GitHub/PhD_IO/pset1/OTC_Headache_upd.csv", header = TRUE)
store_dem <- read.csv("store_demographics.csv")


# Reshaping so that each row corresponds to a week-store-brand-size combination

product_data <- OTC_headache %>%
  pivot_longer(sales_1:prom_11,
    names_to = c(".value", "product-brand"),
    names_pattern = "([A-Za-z]+)[^0-9]+([0-9]+)"
  ) %>%
  left_join(store_dem, by = c("store" = "STORE")) %>%
  clean_names() %>%
  mutate(
    store = as.factor(store),
    product_brand = as.numeric(product_brand)
  ) %>%
  # categorizing products, I assume 1-11 are in the same order as the table in
  # question, and I've bundled all store brand objects as a single brand
  mutate(
    brand_name = case_when(
      product_brand %in% 1:3 ~ "Tylenol",
      product_brand %in% 4:6 ~ "Advil",
      product_brand %in% 7:9 ~ "Bayer",
      product_brand %in% 10:11 ~ "Store Brand"
    ) %>% as.factor(.),
    size = case_when(
      product_brand %in% c(1, 4, 7) ~ 25,
      product_brand %in% c(2, 5, 8, 10) ~ 50,
      product_brand %in% c(3, 6, 9, 11) ~ 100
    )
  ) %>%
  group_by(week) %>%
  mutate(
    haus_iv = (n() * mean(price) - price) / (n() - 1)
  ) %>%
  group_by(week, store) %>%
  # assumption - 2% of drug store customers go to it to by headache medicine
  mutate(
    mkt_share = sales / (count * 0.02)
  ) %>%
  mutate(
    mkt_share_inside = sum(mkt_share),
    mkt_share_outside = 1 - sum(mkt_share),
    prom = ifelse(prom > 0, 1, prom)
    # adjust promotion variable
  ) %>%
  # creating brand-by-store dummy, this will make it easier than doing in estimation itself
  mutate(
    brand_by_store = paste(brand_name, "_", store),
    delta_log_share = log(mkt_share) - log(mkt_share_outside)
  ) %>%
  select(-c(count, sales, educ, hsizeavg, nwhite, mkt_share_inside, mkt_share_outside))

#adding iv vars, even the ones we wont use now (only in blp)
product_data <- product_data %>%
  left_join(OTC_instruments %>%
              mutate(store = factor(store)) %>%
              rename(product_brand = brand) %>%
              select(-cost_)) %>%
  na.omit()
  


product_data_long_blp <- product_data %>%
  dummy_cols(select_columns = "brand_name") %>%
  select(-c(brand_by_store, brand_name)) %>%
  select(price, prom, starts_with('brand_'), cost, haus_iv, avoutprice,
         starts_with('pricestore')
         )

  

write_csv(product_data, file = "data_long.csv")
write_csv(product_data_long_blp, file = "data_BLPlong.csv")




models_noIV <- feols(delta_log_share ~ price + prom | csw0(brand_name, brand_by_store),
  data = product_data
)

models_costIV <- feols(delta_log_share ~ prom | csw0(brand_name, brand_by_store) |
  price ~ cost,
data = product_data
)


models_hausIV <- feols(delta_log_share ~ prom | csw0(brand_name, brand_by_store) |
  price ~ haus_iv,
data = product_data
)


models <- list(
  models_noIV[[1]],
  models_noIV[[2]],
  models_noIV[[3]],
  models_costIV[[1]],
  models_costIV[[2]],
  models_costIV[[3]],
  models_hausIV[[1]],
  models_hausIV[[2]],
  models_hausIV[[3]]
)


x <- modelsummary(models,
  output = "latex",
  statistic = NULL,
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  coef_omit = "Intercept",
  gof_omit = "DF|Deviance|R2|AIC|BIC|RMSE|Std.Errors"
) %>%
  add_header_above(c(" " = 1, "No IV" = 3, "Cost IV" = 3, "Hausman IV" = 3))

save_kable(x, "basic_models.tex")





## Calculating elasticities

# Getting price coefficients from each model


alpha <- models_noIV %>% map(~ .x$coefficients[[2]])



# getting mean price and market-shares (unweighted) by product
mean_price_share <- product_data %>%
  group_by(product_brand) %>%
  summarise(
    mkt_share = mean(mkt_share),
    price = mean(price)
  )

own_price_elasticity <- mean_price_share %>%
  ungroup() %>%
  mutate(
    model_1 = -alpha[[1]] * price * (1 - mkt_share),
    model_2 = -alpha[[2]] * price * (1 - mkt_share),
    model_3 = -alpha[[3]] * price * (1 - mkt_share)
  ) %>%
  left_join(product_data %>%
    ungroup() %>%
    select(product_brand, brand_name, size) %>%
    unique()) %>%
  relocate(brand_name, size) %>%
  select(-c(product_brand, mkt_share, price)) %>%
  rename(
    Brand = brand_name,
    Size = size
  )

y <- datasummary_df(own_price_elasticity, output = "latex") %>%
  add_header_above(c(" " = 2, "Price elasticities" = 3))

save_kable(y, "elasticities_basic.tex")
