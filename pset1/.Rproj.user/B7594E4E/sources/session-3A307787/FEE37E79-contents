### Transforming the cleaned data set back to wide format - this will be useful 
### for BLP

data_long <- read.csv('data_longFormat.csv')

### delete stuff

data_long <- data_long %>%
  mutate(
    branded = ifelse(brand == "Store Brand", 0, 1)
  ) %>%
  select(week, store, product_brand, price, prom, income, mkt_share, branded)


data_clean_wide <- data_long %>%
  pivot_wider(
    names_from = product_brand,
    values_from = price:branded
  ) %>%
  na.omit()

write.csv(data_clean_wide, 'data_clean_wide.csv', row.names = FALSE)
