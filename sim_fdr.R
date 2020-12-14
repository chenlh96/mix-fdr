rm(list = ls())

library(locfdr)
source("mixfdr_fun.R")

set.seed(2020)

## simulate data
n_rep = 100
N = 1000
N_null = 50
J = 3

mu_z = c(rep(0, N - N_null), runif(N_null, 2, 4))
err_mat = matrix(rnorm(N * n_rep), ncol = N, nrow = n_rep)
z_mat = t(mu_z + t(err_mat))
zz = seq(min(z_mat), max(z_mat), by = 0.01)
n_zz = length(zz)

mixfdr_fdrval <- mixfdr_theo_fdrval <-
  locfdr_cme_fdrval <- locfdr_mle_fdrval <- locfdr_theo_fdrval <- matrix(0, nrow = n_rep, ncol = n_zz)

mixfdr_FDRval <- mixfdr_theo_FDRval <-
  locfdr_cme_FDRval <- locfdr_mle_FDRval <- locfdr_theo_FDRval <- matrix(0, nrow = n_rep, ncol = n_zz)

for (i in 1:n_rep) {
  # get data
  z = z_mat[i,]
  
  # mixfdr emp null
  mixfdr_res = mixfdr(z, N, J = J, P = 50, turning_P = FALSE, near_null_thres = 0.1, plot = FALSE)
  mix_prior = mixfdr_res$prior_density
  null_cpns = mixfdr_res$null_cpns
  fdrval = calc_fdr(zz, n_zz, J, mix_prior, null_cpns)
  mixfdr_fdrval[i, ] = fdrval$fdr
  mixfdr_FDRval[i, ] = fdrval$Fdr_2side
  
  # mixfdr theo null
  mixfdr_theo_res = mixfdr(z, N, J = J, P = 50, turning_P = FALSE, theo_null = TRUE, near_null_thres = 0.1, plot = FALSE)
  theo_mix_prior = mixfdr_theo_res$prior_density
  theo_null_cpns = mixfdr_theo_res$null_cpns
  theo_fdrval = calc_fdr(zz, n_zz, J, theo_mix_prior, theo_null_cpns)
  mixfdr_theo_fdrval[i, ] = theo_fdrval$fdr
  mixfdr_theo_FDRval[i, ] = theo_fdrval$Fdr_2side
  
  # locfdr emp null estimated by mle
  locfdr_mle_res = locfdr(z, nulltype = 1, plot = FALSE)
  locfdr_mle_mat = locfdr_mle_res$mat
  locfdr_mle_FDR = cbind(locfdr_mle_mat[,"x"], 
                         apply(locfdr_mle_mat[,c("f", "f0", "f0theo")], 2, function(x) cumsum(x) / sum(x)))
  colnames(locfdr_mle_FDR) <- c("x", "F", "F0", "F0theo")
  F_zz = approx(locfdr_mle_FDR[,'x'], locfdr_mle_FDR[,"F"], zz, rule = 2, ties="ordered")$y
  F_negzz = approx(locfdr_mle_FDR[,'x'], locfdr_mle_FDR[,"F"], -zz, rule = 2, ties="ordered")$y
  F0_mle_zz = approx(locfdr_mle_FDR[,'x'], locfdr_mle_FDR[,"F0"], zz, rule = 2, ties="ordered")$y
  F0_mle_negzz = approx(locfdr_mle_FDR[,'x'], locfdr_mle_FDR[,"F0"], -zz, rule = 2, ties="ordered")$y
  p0_mle = locfdr_mle_res$fp0["mlest", "p0"]
  
  locfdr_mle_fdrval[i, ] = approx(locfdr_mle_mat[,'x'], locfdr_mle_mat[,"fdr"], zz, rule = 2, ties="ordered")$y
  locfdr_mle_FDRval[i, ] = p0_mle * (1 - F0_mle_zz + F0_mle_negzz) / (1 - F_zz + F_negzz)
  locfdr_mle_FDRval[i, ] = pmax(pmin(locfdr_mle_FDRval[i, ], 1), 0)
  
  # locfdr emp null estimated by cme
  locfdr_cme_res = locfdr(z, nulltype = 2, plot = FALSE)
  locfdr_cme_mat = locfdr_cme_res$mat
  locfdr_cme_FDR = cbind(locfdr_cme_mat[,"x"], 
                         apply(locfdr_cme_mat[,c("f", "f0", "f0theo")], 2, function(x) cumsum(x) / sum(x)))
  colnames(locfdr_cme_FDR) <- c("x", "F", "F0", "F0theo")
  F0_cme_zz = approx(locfdr_cme_FDR[,'x'], locfdr_cme_FDR[,"F0"], zz, rule = 2, ties="ordered")$y
  F0_cme_negzz = approx(locfdr_cme_FDR[,'x'], locfdr_cme_FDR[,"F0"], -zz, rule = 2, ties="ordered")$y
  p0_cme = locfdr_cme_res$fp0["mlest", "p0"]
  
  locfdr_cme_fdrval[i, ] = approx(locfdr_cme_mat[,'x'], locfdr_cme_mat[,"fdrtheo"], zz, rule = 2, ties="ordered")$y
  locfdr_cme_FDRval[i, ] = p0_cme * (1 - F0_cme_zz + F0_cme_negzz) / (1 - F_zz + F_negzz)
  locfdr_cme_FDRval[i, ] = pmax(pmin(locfdr_cme_FDRval[i, ], 1), 0)
  
  # locfdr theo null
  # already been computed/existed in each locfdr output list
  F0_theo_zz = approx(locfdr_mle_FDR[,'x'], locfdr_mle_FDR[,"F0theo"], zz, rule = 2, ties="ordered")$y
  F0_theo_negzz = approx(locfdr_mle_FDR[,'x'], locfdr_mle_FDR[,"F0theo"], -zz, rule = 2, ties="ordered")$y
  p0_theo = locfdr_cme_res$fp0["thest", "p0"]
  
  locfdr_theo_fdrval[i, ] = approx(locfdr_mle_mat[,'x'], locfdr_mle_mat[,"fdrtheo"], zz, rule = 2, ties="ordered")$y
  locfdr_theo_FDRval[i, ] = p0_theo * (1 - F0_theo_zz + F0_theo_negzz) / (1 - F_zz + F_negzz)
  locfdr_theo_FDRval[i, ] = pmax(pmin(locfdr_theo_FDRval[i, ], 1), 0)
}

J_true = 2
true_params = cbind("piprob" = c(950, 50) / 1000,
                    "mu" = c(0, mean(mu_z[951:1000])),
                    "sigma" = c(0.001, var(mu_z[951:1000])))
truefdr_res = calc_fdr(zz, n_zz, J_true, true_params, c(TRUE, FALSE))
truefdr = truefdr_res$fdr
trueFDR = truefdr_res$Fdr_2side

mixfdr_fdr_mean = colMeans(mixfdr_fdrval)
mixfdr_fdr_sd = apply(mixfdr_fdrval, 2, sd)
mixfdr_FDR_mean = colMeans(mixfdr_FDRval)
mixfdr_FDR_sd = apply(mixfdr_FDRval, 2, sd)

mixfdr_theo_fdr_mean = colMeans(mixfdr_theo_fdrval)
mixfdr_theo_fdr_sd = apply(mixfdr_theo_fdrval, 2, sd)
mixfdr_theo_FDR_mean = colMeans(mixfdr_theo_FDRval)
mixfdr_theo_FDR_sd = apply(mixfdr_theo_FDRval, 2, sd)

locfdr_theo_fdr_mean = colMeans(locfdr_theo_fdrval)
locfdr_theo_fdr_sd = apply(locfdr_theo_fdrval, 2, sd)
locfdr_theo_FDR_mean = colMeans(locfdr_theo_FDRval)
locfdr_theo_FDR_sd = apply(locfdr_theo_FDRval, 2, sd)

locfdr_mle_fdr_mean = colMeans(locfdr_mle_fdrval)
locfdr_mle_fdr_sd = apply(locfdr_mle_fdrval, 2, sd)
locfdr_mle_FDR_mean = colMeans(locfdr_mle_FDRval)
locfdr_mle_FDR_sd = apply(locfdr_mle_FDRval, 2, sd)

locfdr_cme_fdr_mean = colMeans(locfdr_cme_fdrval)
locfdr_cme_fdr_sd = apply(locfdr_cme_fdrval, 2, sd)
locfdr_cme_FDR_mean = colMeans(locfdr_cme_FDRval)
locfdr_cme_FDR_sd = apply(locfdr_cme_FDRval, 2, sd)

cols = rainbow(5)

par(mfrow = c(1, 2))
plot(zz, truefdr, 'l', col = "black", lty = 1, lwd = 2,
     xlab = 'z', ylab = "fdr", main = "Efdr")
lines(zz, mixfdr_fdr_mean, lty = 1, col = cols[1])
lines(zz, mixfdr_theo_fdr_mean, lty = 1, col = cols[2])
lines(zz, locfdr_theo_fdr_mean, lty = 2, col = cols[3])
lines(zz, locfdr_cme_fdr_mean, lty = 2, col = cols[4])
lines(zz, locfdr_mle_fdr_mean, lty = 2, col = cols[5])
legend("topright", legend = c("true", "mix emp", "mix theo", "loc theo", "loc cme", "loc mle"), 
       lty = c(1, 1,1,3,3,3), col = c("black", cols), bty = 'n', cex = 0.6)

plot(zz, mixfdr_fdr_sd, 'l', col = cols[1],
     xlab = 'z', ylab = "fdr", main = "Var fdr")
lines(zz, mixfdr_theo_fdr_sd, lty = 1, col = cols[2])
lines(zz, locfdr_theo_fdr_sd, lty = 3, col = cols[3])
lines(zz, locfdr_cme_fdr_sd, lty = 3, col = cols[4])
lines(zz, locfdr_mle_fdr_sd, lty = 3, col = cols[5])
legend("topright", legend = c("mix emp", "mix theo", "loc theo", "loc cme", "loc mle"), lty = c(1,1,3,3,3), col = cols,
       bty = 'n', cex = 0.6)

par(mfrow = c(1, 2))
plot(zz, trueFDR, 'l', col = "black", lty = 1, lwd = 2,
     xlab = 'z', ylab = "fdr", main = "Efdr")
lines(zz, mixfdr_FDR_mean, lty = 1, col = cols[1])
lines(zz, mixfdr_theo_FDR_mean, lty = 1, col = cols[2])
lines(zz, locfdr_theo_FDR_mean, lty = 2, col = cols[3])
lines(zz, locfdr_cme_FDR_mean, lty = 2, col = cols[4])
lines(zz, locfdr_mle_FDR_mean, lty = 2, col = cols[5])
legend("topright", legend = c("true", "mix emp", "mix theo", "loc theo", "loc cme", "loc mle"), 
       lty = c(1, 1,1,3,3,3), col = c("black", cols), bty = 'n', cex = 0.6)

plot(zz, mixfdr_FDR_sd, 'l', col = cols[1],
     xlab = 'z', ylab = "fdr", main = "Var fdr")
lines(zz, mixfdr_theo_FDR_sd, lty = 1, col = cols[2])
lines(zz, locfdr_theo_FDR_sd, lty = 3, col = cols[3])
lines(zz, locfdr_cme_FDR_sd, lty = 3, col = cols[4])
lines(zz, locfdr_mle_FDR_sd, lty = 3, col = cols[5])
legend("topright", legend = c("mix emp", "mix theo", "loc theo", "loc cme", "loc mle"), lty = c(1,1,3,3,3), col = cols,
       bty = 'n', cex = 0.6)


