rm(list = ls())

library(splines)
library(locfdr)
source("mixfdr_fun.R")

set.seed(10)

n_rep = 100
mu_all = 2:5
K_all = 5 * c(1, 10, 100)
N = 1000

params_all = expand.grid(mu = mu_all, K = K_all)
gen_data_all = vector("list", nrow(params_all))

effect_size_spline <- function(z,df=5){
    bins = seq(min(z)-.1,max(z)+.1, len = 120)
    h = hist(z, bins, plot = F)
    x = h$mids
    y = h$counts
    g  = glm(y ~ ns(x,df=df),family = "poisson")
    ss = splinefun(x, log(dnorm(x) / g$fit),method = "natural")
    return(-ss(z,deriv=1))		
  }

# spline_df5_mse0 <- bayes_mse0 <- ebthres_mse0 <- mixfdr_mse0 <- rep(0, nrow(params_all))
# for (p in 1:nrow(params_all)) {
#   mu = params_all$mu[p]
#   K = params_all$K[p]
#   mu_z = c(rep(0, N - K), runif(K, mu - 0.5, mu + 0.5))
#   err_mat = matrix(rnorm(N * n_rep), nrow = n_rep, ncol = N)
#   z_mat = t(mu_z + t(err_mat))
#   J_true = 2
#   true_prior = cbind("piporb" = c(N - K, K) / N, 
#                      "mu" = c(0, mu),
#                      "sigma2" = c(0.001, 1 / 12))
#   
#   ebthres_est <- spline_est <-
#     bayes_est <- mixfdr_est <- matrix(0, nrow = n_rep, ncol = N)
#   
#   for (i in 1:n_rep) {
#     z = z_mat[i, ]
#     
#     mixfdr_res = mixfdr(z, N, J = 9, P = 50, turning_P = FALSE, near_null_thres = 0.1, plot = FALSE)
#     mixfdr_est[i, ] = mixfdr_res$effect_sz
#     bayes_est[i, ] = calc_effect_size(z, N, J_true, true_prior)$effect_sz
#     ebthres_est[i, ] = ebayesthresh(z)
#     spline_est[i, ] = effect_size_spline(z)
#   }
#   
#   mixfdr_mse0[p] = mean(rowSums(t(mu_z - t(mixfdr_est))^2))
#   bayes_mse0[p] = mean(rowSums(t(mu_z - t(bayes_est))^2))
#   ebthres_mse0[p] = mean(rowSums(t(mu_z - t(ebthres_est))^2))
#   spline_df5_mse0[p] = mean(rowSums(t(mu_z - t(spline_est))^2))
# }
# 
# mixfdr_mse = matrix(mixfdr_mse0, nrow = length(K_all), ncol = length(mu_all), byrow = TRUE)
# bayes_mse = matrix(bayes_mse0, nrow = length(K_all), ncol = length(mu_all), byrow = TRUE)
# ebthres_mse = matrix(ebthres_mse0, nrow = length(K_all), ncol = length(mu_all), byrow = TRUE)
# spline_df5_mse = matrix(spline_df5_mse0, nrow = length(K_all), ncol = length(mu_all), byrow = TRUE)
# par(mfrow = c(1, length(K_all)))
# for (i in 1:length(K_all)) {
#   plot(mu_all, mixfdr_mse[i, ] / bayes_mse[i, ], type = "l", lty = 1, col = "black", ylim = c(1, 6),
#        xlab = "mu", ylab = "mse", main = sprintf("K = %d", K_all[i]))
#   lines(mu_all, bayes_mse[i, ] / bayes_mse[i, ], type = "l", lty = 2, col = "red")
#   lines(mu_all, ebthres_mse[i, ] / bayes_mse[i, ], type = "l", lty = 2, col = "blue")
#   lines(mu_all, spline_df5_mse[i, ] / bayes_mse[i, ], type = "l", lty = 2, col = "green")
#   legend("topleft", legend = c("mixfdr", "bayes", "EbayesThresh", "spline"), lty = c(1,2,2,2), 
#          col = c("black", "red", "blue", "green"), bty = "n", cex = 1)
# }

spline_df5_mse1 <- bayes_mse1 <- ebthres_mse1 <- mixfdr_mse1 <- rep(0, nrow(params_all))
for (p in 1:nrow(params_all)) {
  mu = params_all$mu[p]
  K = params_all$K[p]
  mu_z = c(rep(0, N - K), runif(ceiling(K * 2 / 3), mu - 0.5, mu + 0.5), runif(K - ceiling(K * 2 / 3), -mu - 0.5, -mu + 0.5))
  err_mat = matrix(rnorm(N * n_rep), nrow = n_rep, ncol = N)
  z_mat = t(mu_z + t(err_mat))
  J_true = 3
  true_prior = cbind("piporb" = c(N - K, ceiling(K * 2 / 3), K - ceiling(K * 2 / 3)) / N, 
                     "mu" = c(0, mu, -mu),
                     "sigma2" = c(0.001, 1 / 12, 1 / 12))
  
  ebthres_est <- spline_est <-
    bayes_est <- mixfdr_est <- matrix(0, nrow = n_rep, ncol = N)
  
  for (i in 1:n_rep) {
    z = z_mat[i, ]
    
    mixfdr_res = mixfdr(z, N, J = 9, P = 50, turning_P = FALSE, near_null_thres = 0.1, plot = FALSE)
    mixfdr_est[i, ] = mixfdr_res$effect_sz
    bayes_est[i, ] = calc_effect_size(z, N, J_true, true_prior)$effect_sz
    ebthres_est[i, ] = ebayesthresh(z)
    spline_est[i, ] = effect_size_spline(z)
  }
  
  mixfdr_mse1[p] = mean(rowSums(t(mu_z - t(mixfdr_est))^2))
  bayes_mse1[p] = mean(rowSums(t(mu_z - t(bayes_est))^2))
  ebthres_mse1[p] = mean(rowSums(t(mu_z - t(ebthres_est))^2))
  spline_df5_mse1[p] = mean(rowSums(t(mu_z - t(spline_est))^2))
}

mixfdr_mse = matrix(mixfdr_mse1, nrow = length(K_all), ncol = length(mu_all), byrow = TRUE)
bayes_mse = matrix(bayes_mse1, nrow = length(K_all), ncol = length(mu_all), byrow = TRUE)
ebthres_mse = matrix(ebthres_mse1, nrow = length(K_all), ncol = length(mu_all), byrow = TRUE)
spline_df5_mse = matrix(spline_df5_mse1, nrow = length(K_all), ncol = length(mu_all), byrow = TRUE)
par(mfrow = c(1, length(K_all)))
for (i in 1:length(K_all)) {
  plot(mu_all, mixfdr_mse[i, ] / bayes_mse[i, ], type = "l", lty = 1, col = "black", ylim = c(1, 6),
       xlab = "mu", ylab = "mse", main = sprintf("K = %d", K_all[i]))
  lines(mu_all, bayes_mse[i, ] / bayes_mse[i, ], type = "l", lty = 2, col = "red")
  lines(mu_all, ebthres_mse[i, ] / bayes_mse[i, ], type = "l", lty = 2, col = "blue")
  lines(mu_all, spline_df5_mse[i, ] / bayes_mse[i, ], type = "l", lty = 2, col = "green")
  legend("topleft", legend = c("mixfdr", "bayes", "EbayesThresh", "spline"), lty = c(1,2,2,2), 
         col = c("black", "red", "blue", "green"), bty = "n", cex = 1)
}

