
library(tidyverse)
library(gmm)


set.seed(2)

T_ <-  10
I <-  1000


# simulate productivity 50 periods past just to get stationarity
omega_past <-  matrix(0, I, 50)

for (t in 2:50){
  omega_past[, t] = 0.9*omega_past[, t-1] + rnorm(I, 0, 1)
}

xi <-  matrix(rnorm(I*T_, 0, 1), I, T_)

omega <-  matrix(0, I, T_)
omega[, 1] = 0.9*omega_past[, 50] + xi[, 1]

for (t in 2:T_){
  omega[, t] = 0.9*omega[, t-1] + xi[, t]
}

epsilon <-  matrix(rnorm(I*T_, 0, 3), I, T_)
u <-  matrix(rnorm(I*T_, 0, 0.2), I, T_)
# capital and labor are determined one period in advance

L <-  matrix(0, I, T_)
L[, 1] <-  (0.108)^3*(0.73)^(0.8)*exp(3.69 + 3.6*omega_past[,50])
L[, 2:10] <-  (0.108)^3*(0.73)^(0.8)*exp(3.69 + 3.6*omega[,1:9])

K = 0.73*L

# materials are determined on the same period, and can be conditioned on epsilon

M <-  (0.4)^(5/3)*exp(8/75)*(exp(0.3*epsilon + omega)*L^(0.3)*K^(0.2))^(4/3)

# Realized revenue depends on the actual realizations of u

R <-  (exp(u + omega + 0.3*epsilon)*L^(0.3)*K^(0.2)*M^(0.5))^(0.8) 


rm(epsilon, omega, omega_past, u, xi, t)
# Putting in long format


K = data.frame(K) 
L = data.frame(L)
M = data.frame(M)
R = data.frame(R)

l1 = list(K = K, L = L, M = M, R = R)
l2 = list('K', 'L', 'M', 'R')

dataset <- 
  map2(l1, l2, ~ .x %>%
        mutate(i = row_number()) %>%
        relocate(i) %>%
        pivot_longer(X1:X10, names_to = "t",
                     values_to = .y) %>%
         mutate(t = substr(t, 2, nchar(t)) %>%
                  as.numeric()))


data_long <-  cbind(dataset$K, 
                  dataset$L[, 3],
                  dataset$M[,3],
                  dataset$R[, 3])

# getting lagged inputs
data_long <- data_long %>%
  group_by(i) %>%
  arrange(t) %>%
  mutate(L_lag = lag(L, 1),
         K_lag = lag(K, 1),
         M_lag = lag(M, 1)) %>%
  na.omit()


write_csv()





