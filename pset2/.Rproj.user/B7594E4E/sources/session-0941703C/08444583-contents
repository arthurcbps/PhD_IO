
library(tidyverse)

library(gmm)


moment <- function(par, x){
  
  beta_0 = par[1]
  beta_k = par[2]  
  beta_l = par[3]
  rho = par[4]
  
  k = log(x[, 1])
  l = log(x[, 2]) 
  m = log(x[, 3])
  r = log(x[, 4])
  lag_k = log(x[, 5])
  lag_l = log(x[, 6])
  lag_m = log(x[, 7])
  phi_hat = x[, 8]
  lag_phi_hat = x[, 9]
  
  xi = phi_hat - beta_0 - beta_k*k - beta_l*l - rho*(lag_phi_hat - beta_0 - beta_k*lag_k - beta_l*lag_l)
  
  m1 = xi
  m2 = xi*k
  m3 = xi*l
  m4 = xi*lag_phi_hat
  
  f = cbind(m1, m2, m3, m4)
  
  return(f)
  
}

gmm_crit <- function(par){
  x <- moment(par, x = data_long %>%
                ungroup() %>%
                select(K, L, M, R, K_lag, L_lag, M_lag, phi_hat, phi_hat_lag) %>%
                as.matrix())
  
  
  y <- matrix(colMeans(x)) 
  yt <- t(colMeans(x))
  
  return(yt%*%y)
}


T_ <-  10
I <-  1000
D = 1


rho_hat <- rep(0, D)



for (d in 1:D) {
  

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

rm(l1, l2, K, L, M, R)

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
         M_lag = lag(M, 1),
         R_lag = lag(R, 1)) 


#write_csv(data_long, file = "sim_dgp_data.csv")


# ACF 2 step procedure

# 1st step is the filtering stage a nonparametric regression of log revenue on inputs - 
# we use then a 4th degree polynomial


first_stage <- lm(log(R) ~ polym(log(K), log(L), log(M), degree = 2),
                  data = data_long)


data_long$phi_hat = fitted(first_stage)

data_long <- data_long %>%
  group_by(i) %>%
  arrange(t) %>%
  mutate(phi_hat_lag = lag(phi_hat, 1)) %>%
  na.omit()

# second stage - GMM




acf <- optim(gmm_crit,
             par = c(beta_0 = 0.1, beta_k = 1, beta_l = 1.5, rho = 0.9))

rho_hat[d] <- acf$par[4]

}



