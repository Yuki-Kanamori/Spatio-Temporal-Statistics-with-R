#------------------------------------------
# Supplyment of the first part of Chap.4
# developed by Tetsuya Akita
#
# 1. Estimating Y and sigma with known beta
# 2. Estimating Y and sigma with unknown beta
# 3. Replicate Fig.4.5
#------------------------------------------

library(tidyverse)

# 1. Estimating Y and sigma with known beta

Z = tibble::tibble(s=c(2,2,6,6),t=c(0.2,1,0.2,0.9),n=c(15,22,17,23))
Y = tibble::tibble(s=c(3),t=c(0.5),n=c(NA))
cov_func = function(a=2,b=0.2,v=2,d=1,s,t){
  v * exp(-b^2*s^2/(a^2*t^2+1)) / (a^2*t^2+1)^(d/2)
}

Cz = as.matrix(cov_func(s=dist(Z$s,diag = T,upper = T),
                        t=dist(Z$t,diag = T,upper = T))) + 2*diag(4) # cov_func does not return the diag elements...
c0 = cov_func(s=Z$s-Y$s, t=Z$t-Y$t)
w_t = t(c0) %*% solve(Cz)

beta = 20
X = c(1,1,1,1)
cond_mean = as.vector(t(Z$n) - X*beta)
Y_hat = beta + t(c0) %*% solve(Cz) %*% cond_mean
sigma_sk = 2- t(c0) %*% solve(Cz) %*% c0
cat("Y_hat =",Y_hat,"sigma_sk =",sigma_sk,"\n")

ggplot(data=data.frame(X_axis=c(15,25)), aes(x=X_axis)) +
  stat_function(fun=dnorm, args=list(mean=beta, sd=2^0.5)) +
  stat_function(fun=dnorm, args=list(mean=Y_hat, sd=sigma_sk),col="red")

# 2. Estimating Y and sigma with unknown beta

beta_hat = solve( t(X) %*% solve(Cz) %*% X ) %*% t(X) %*% solve(Cz) %*% Z$n
cat("beta_hat =",beta_hat,"\n")

X = c(1,1,1,1)
cond_mean = as.vector(t(Z$n) - X*c(beta_hat))
Y_hat_uk = beta_hat + t(c0) %*% solve(Cz) %*% cond_mean
kap = t(1 - t(X) %*% solve(Cz) %*% c0) %*% solve( t(X) %*% solve(Cz) %*% X ) %*% (1 - t(X) %*% solve(Cz) %*% c0)
sigma_uk = 2- t(c0) %*% solve(Cz) %*% c0 + kap
cat("Y_hat_uk =",Y_hat_uk,"sigma_uk =",sigma_uk,"\n")

ggplot(data=data.frame(X_axis=c(15,25)), aes(x=X_axis)) +
  stat_function(fun=dnorm, args=list(mean=beta, sd=2^0.5)) +
  stat_function(fun=dnorm, args=list(mean=Y_hat, sd=sigma_sk),col="red") +
  stat_function(fun=dnorm, args=list(mean=Y_hat_uk, sd=sigma_uk),col="blue")

# 3. Replicate Fig.4.5
rm(list = ls())
library(tidyverse)

assign_condtion = function(n_max,t_max){
  out = NULL
  ran = runif(n_max)
  for(i in 1:n_max){
    if (ran[i]<1/3) {
      out = c(out,rep("control",t_max))
    } else if (ran[i] < 2/3) {
      out = c(out,rep("treatment_1",t_max))
    } else {
      out = c(out,rep("treatment_2",t_max))
    }
  }
  return(out)
}

n_max = 90
t_max = 20
sigma_1 = 0.2
sigma_2 = 0.2
sigma_ep = 0.1
b_0 = 25
b_1 = 2
b_2 = 2.5
b_3 = 3

param = tibble::tibble(
  n = seq(1,n_max),
  alpha_1 = rnorm(n = n_max, mean = 0, sd = sigma_1^0.5),
  alpha_2 = rnorm(n = n_max, mean = 0, sd = sigma_2^0.5),
  intecept =  b_0 + rnorm(n = n_max, mean = 0, sd = sigma_1^0.5)
) 
print(param)

df = tidyr::crossing(
  n = seq(1,n_max), 
  t = seq(1,t_max)
) %>%
  dplyr::mutate(condition = assign_condtion(n_max,t_max)) %>%
  dplyr::mutate(epsilon =  rnorm(n = n_max*t_max, mean = 0, sd = sigma_ep^0.5)) %>%
  dplyr::mutate(intercept = purrr::map_dbl(n, ~ {param$intecept[.]})) %>%
  dplyr::mutate(slope = if_else(condition == "control",
                                true = purrr::map2_dbl(n, t, function(x,y){
                                  (b_1 + param$alpha_2[x]) * y}),
                                false = if_else(condition == "treatment_1",
                                                true = purrr::map2_dbl(n, t, function(x,y){
                                                  (b_2 + param$alpha_2[x]) * y
                                                  }),
                                                false = purrr::map2_dbl(n, t, function(x,y){
                                                  (b_3 + param$alpha_2[x]) * y
                                                })))) %>%
  dplyr::mutate(Z = intercept + slope + epsilon)
print(df)
ggplot(df,aes(x=t,y=Z,colour=condition,group=n)) +
  geom_line()







