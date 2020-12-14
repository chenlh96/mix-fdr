gen_random_data <- function(N, J, piprob, mu, sigma2, only_obs = TRUE) {
  id_cpn = apply(rmultinom(N, 1, piprob), 2, which.max)
  mu_z = rnorm(N, mean = mu[id_cpn], sd = sqrt(sigma2[id_cpn]))
  z = rnorm(N, mean = mu_z)
  if (only_obs) return(z) ## dirty implementation, violate the principle that
                          ## output should retain the same class
  dat = list(mu_z = mu_z, z = z, params = cbind("piprob" = piprob, "mu" = mu, "sigma2" = sigma2))
  return(dat)
}

turning_penalty_param <- function(z, N, J, K, L, B, theo_null = FALSE, keep_gen_data = FALSE) {
  pen_candid = ceiling(seq(100, N / 2, length.out = K))
  pen_ori = N / 5
  
  # pre_model = mix_prior_EM(z, N, J, pen_ori, theo_null)
  pre_model = mix_prior_EM_multi_start(z, N, J, pen_ori, theo_null)
  nz = matrix(rnorm(J* L * 2, sd = 0.1), nrow = L)
  pert_model = vector("list", L)
  for (l in 1:L) {
    nz_mu_sig = matrix(c(rep(0, J), nz[l,]), nrow = J)
    pert_params = pre_model$params + nz_mu_sig
    pert_params[, 3] = ifelse(pert_params[, 3] <= 0, rep(1, J), pert_params[, 3])
    pert_model[[l]] = list(N = N, J = J, params = pert_params)
  }
  
  rand_dataset = matrix(0, nrow = N, ncol = B * L)
  rand_pert_params = vector("list", B * L)
  for (b in 1:B) {
    for (l in 1:L) {
      rand_dataset[,(b - 1) * 10 + l] = 
        gen_random_data(N, J, pert_model[[l]]$params[,1], pert_model[[l]]$params[,2], pert_model[[l]]$params[,3])
      rand_pert_params[[(b - 1) * 10 + l]] = pert_model[[l]]$params
    }
  }
  
  cv_loss <- matrix(0, ncol = K, nrow = B * L)
  for (i in 1:(B * L)) {
    for (k in 1:K) {
      pert_fit = mix_prior_EM_multi_start(rand_dataset[,i], N, J, pen_candid[k], theo_null)
      cv_loss[i, k] = sum((pert_fit$params - rand_pert_params[[i]])^2)
    }
  }
  pen_best_idx = which.max(colMeans(cv_loss))
  
  P = pen_candid[pen_best_idx]
  res = list(P = P, cv_loss = cv_loss, P_candidate = pen_candid, P_ori = pen_ori,
             pre_model = pre_model, pert_model = pert_model)
  if (keep_gen_data) res[["random_dataset"]] = rand_dataset
  return(res)
}

mix_prior_EM_multi_start <- function(z, N, J, P, theo_null = FALSE, max_iter = 200, tol = 1e-3) {
  # first do default start
  default_res = mix_prior_EM(z, N, J, P, theo_null)
  
  n_aux = 3
  center_area = matrix(c(0.2, 0.8,
                         0.15, 0.85,
                         0.1, 0.9), nrow =n_aux, ncol = 2, byrow = TRUE)
  aux_start = array(0, c(3, J, 3))
  for (i in 1:n_aux) {
    split_cpn0 = center_area[i,]
    prop_non0 = (1 - diff(split_cpn0)) / (J - 1)
    split_cpn = ceiling(c(seq(0, split_cpn0[1], by = prop_non0), seq(split_cpn0[2], 1, by = prop_non0)) * N)
    split_idx = numeric(N)
    for(j in 1:J) split_idx[(split_cpn[j] + 1):split_cpn[j + 1]] = j
    z_group = split(sort(z), split_idx)
    j0 = ceiling(J / 2)
    z_group = c(z_group[j0], z_group[-j0])
    names(z_group) = 1:J
    
    piprob0 = vapply(z_group, function(x) length(x) / N, numeric(1))
    mu0 = vapply(z_group, mean, numeric(1))
    sigma20 = vapply(z_group, function(x) max(c(var(x) - 1), 0), numeric(1))
    aux_start[i,,] = cbind(piprob0, mu0, sigma20)
  }
  
  multi_start_res = vector("list", n_aux)
  loglik = rep(0, n_aux)
  for (i in 1:n_aux) {
    multi_start_res[[i]] = mix_prior_EM(z, N, J, P, theo_null, aux_start[i,,])
    loglik[i] = multi_start_res[[i]]$loglik[multi_start_res[[i]]$stop_iter - 1]
  }
  
  loglik = c(default_res$loglik[default_res$stop_iter - 1], loglik)
  best_idx = which.max(loglik)
  if (best_idx == 1) {
    res = default_res
  } else res = multi_start_res[[best_idx - 1]]
  
  return(res)
}

mix_prior_EM <- function(z, N, J, P, theo_null = FALSE, init_params = NULL, max_iter = 200, tol = 1e-3) {
  ## set starting point
  if (is.null(init_params)) {
    ## initializing start point
    split_cpn0 = c(0.25, 0.75)
    prop_non0 = (1 - diff(split_cpn0)) / (J - 1)
    split_cpn = ceiling(c(seq(0, split_cpn0[1], by = prop_non0), seq(split_cpn0[2], 1, by = prop_non0)) * N)
    split_idx = numeric(N)
    for(i in 1:J) split_idx[(split_cpn[i] + 1):split_cpn[i + 1]] = i
    z_group = split(sort(z), split_idx)
    j0 = ceiling(J / 2)
    z_group = c(z_group[j0], z_group[-j0])
    names(z_group) = 1:J
    
    piprob <- piprob0 <- vapply(z_group, function(x) length(x) / N, numeric(1))
    mu <- mu0 <- vapply(z_group, mean, numeric(1))
    sigma2 <- sigma20 <- vapply(z_group, function(x) max(c(var(x) - 1), 0), numeric(1))
  } else {
    if (!("matrix" %in% class(init_params))) stop("not a matrix")
    if (NROW(init_params) != J) stop("not matching number of components")
    if (NCOL(init_params) != 3) stop("not matching number of parameters per component")
    
    piprob <- piprob0 <- init_params[,1]
    mu <- mu0 <- init_params[,2]
    sigma2 <- sigma20 <- init_params[,3]
  }
  pen_beta = c(P, rep(0, J- 1))
  
  
  if (theo_null) {
    mu[1] <- mu0[1] <- 0
    sigma2[1] <- sigma20[1] <- 0
  }
  
  ## EM part
  z_cpn_prob <- pz <- matrix(0, nrow = N, ncol = J)
  loglik = numeric(max_iter)
  loglik[1] = -Inf
  for (i in 2:max_iter) {
    # E-step
    for (j in 1:J) z_cpn_prob[,j] = piprob[j] * dnorm(z, mean =mu[j], sd = sqrt(sigma2[j] + 1))
    pz = z_cpn_prob / rowSums(z_cpn_prob)
    
    # check margin lik
    pen_loglik = lgamma(sum(pen_beta[1])) - sum(lgamma(pen_beta[1])) + sum((pen_beta[1] - 1) * log(piprob[1]))
    loglik[i] = sum(log(rowSums(z_cpn_prob))) + pen_loglik
    t = loglik[i] - loglik[i - 1]
    if (is.na(t) | is.nan(t)) t = Inf
    if (t < 0) warning("loglik decrease") #stop("loglik decrease")
    if (t < tol) break
    
    # M-step
    piprob = (colSums(pz) + pen_beta) / (sum(pz) + sum(pen_beta))
    mu = colSums(pz * z) / colSums(pz)
    sigma2 = pmax(colSums(pz * z^2) / colSums(pz) - mu^2 - 1, 0)
    
    if (theo_null) {
      mu[1] = 0
      sigma2[1] = 0
    }
  }
  stop_iter = i
  params = cbind("piprob" = piprob, "mu" = mu, "sigma2" = sigma2)
  params0 = cbind("piprob" = piprob0, "mu" = mu0, "sigma2" = sigma20)
  
  # plot(loglik[2:stop_iter])
  # print(diff(loglik))
  
  res = list(params=params, params0=params0,
             loglik = loglik[2:stop_iter], stop_iter = stop_iter)
  return(res)
}

calc_fdr <- function(z, N, J, mix_prior_den, null_cpns) {
  piprob = mix_prior_den[,1]
  mu = mix_prior_den[,2]
  sigma2 = mix_prior_den[,3]
  
  # estimate fdr, effect size
  z_cpn_prob = matrix(0, nrow = N, ncol = J)
  for (j in 1:J) z_cpn_prob[,j] = piprob[j] * dnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1))
  pz = z_cpn_prob / rowSums(z_cpn_prob)
  prec2 = 1 / sigma2
  z_cpn_psig2 = 1 / (1 + 1/ sigma2)
  z_cpn_pmu = (1 / (1 + prec2)) %o% z + (1 - 1 / (1 + prec2)) * mu
  z_cpn_pmu = t(z_cpn_pmu)
  lfdr = rowSums(pz[,null_cpns,drop = FALSE])
  
  ## estimate the FDR (global fdr,2 side, left and right)
  z_cpn_ldist <- z_cpn_rdist <- z_cpn_twosize <- matrix(0, nrow = N, ncol = J)
  for (j in 1:J) {
    z_cpn_rdist[,j] = piprob[j] * (1 - pnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1)))
    z_cpn_ldist[,j] = piprob[j] * pnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1))
    z_cpn_twosize[,j] = piprob[j] * pnorm(-z, mean = mu[j], sd = sqrt(sigma2[j] + 1)) + 
      piprob[j] * (1 - pnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1)))
  }
  
  FDR = rowSums(z_cpn_twosize[,null_cpns,drop = FALSE]) / rowSums(z_cpn_twosize)
  FDR_left = rowSums(z_cpn_ldist[,null_cpns,drop = FALSE]) / rowSums(z_cpn_ldist)
  FDR_right = rowSums(z_cpn_rdist[,null_cpns,drop = FALSE]) / rowSums(z_cpn_rdist)
  fdr_res = list(fdr = lfdr, Fdr_2side = FDR, Fdr_left = FDR_left, Fdr_right = FDR_right)
  return(fdr_res)
}

calc_effect_size <- function(z, N, J, mix_prior_den) {
  piprob = mix_prior_den[,1]
  mu = mix_prior_den[,2]
  sigma2 = mix_prior_den[,3]
  
  # estimate fdr, effect size
  z_cpn_prob = matrix(0, nrow = N, ncol = J)
  for (j in 1:J) z_cpn_prob[,j] = piprob[j] * dnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1))
  pz = z_cpn_prob / rowSums(z_cpn_prob)
  prec2 = 1 / sigma2
  z_cpn_psig2 = 1 / (1 + 1/ sigma2)
  z_cpn_pmu = (1 / (1 + prec2)) %o% z + (1 - 1 / (1 + prec2)) * mu
  z_cpn_pmu2 = z_cpn_pmu^2 + z_cpn_psig2
  z_cpn_pmu = t(z_cpn_pmu)
  z_cpn_pmu2= t(z_cpn_pmu2)
  
  effect_sz = rowSums(pz * z_cpn_pmu)
  var_effect_sz = rowSums(pz * z_cpn_pmu2) - effect_sz^2
  effsz_res = list(effect_sz = effect_sz, var_effect_sz = var_effect_sz)
  
  return(effsz_res)
}


mixfdr <- function(z, N, J, P = NA, turning_P = TRUE, turning_control = list(K = 10, B = 10, L = 10), theo_null = FALSE, near_null_thres = 0, max_iter = 100, tol = 1e-3, plot = TRUE) {
  if (turning_P) {
    K = turning_control$K
    L = turning_control$L
    B = turning_control$B
    turn_res = turning_penalty_param(z, N, J, K, L, B)
    P = turn_res$P
  } 
  
  mix_prior = mix_prior_EM_multi_start(z, N, J, P, theo_null, max_iter, tol)
  mix_prior_den = mix_prior$params
  piprob = mix_prior_den[,1]
  mu = mix_prior_den[,2]
  sigma2 = mix_prior_den[,3]
  
  # estimate fdr, effect size
  z_cpn_prob = matrix(0, nrow = N, ncol = J)
  for (j in 1:J) z_cpn_prob[,j] = piprob[j] * dnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1))
  pz = z_cpn_prob / rowSums(z_cpn_prob)
  prec2 = 1 / sigma2
  z_cpn_psig2 = 1 / (1 + 1/ sigma2)
  z_cpn_pmu = (1 / (1 + prec2)) %o% z + (1 - 1 / (1 + prec2)) * mu
  z_cpn_pmu2 = z_cpn_pmu^2 + z_cpn_psig2
  z_cpn_pmu = t(z_cpn_pmu)
  z_cpn_pmu2= t(z_cpn_pmu2)
  
  if (!theo_null) null_cpns = abs(mu[1] - mu) <= near_null_thres
  if (theo_null) null_cpns = c(TRUE, rep(FALSE, J - 1))
  lfdr = rowSums(pz[,null_cpns,drop = FALSE])
  effect_sz = rowSums(pz * z_cpn_pmu)
  var_effect_sz = rowSums(pz * z_cpn_pmu2) - effect_sz^2
  
  ## estimate the FDR (global fdr,2 side, left and right)
  z_cpn_ldist <- z_cpn_rdist <- z_cpn_twosize <- matrix(0, nrow = N, ncol = J)
  for (j in 1:J) {
    z_cpn_rdist[,j] = piprob[j] * (1 - pnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1)))
    z_cpn_ldist[,j] = piprob[j] * pnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1))
    z_cpn_twosize[,j] = piprob[j] * pnorm(-z, mean = mu[j], sd = sqrt(sigma2[j] + 1)) + 
      piprob[j] * (1 - pnorm(z, mean = mu[j], sd = sqrt(sigma2[j] + 1)))
  }
  
  FDR = rowSums(z_cpn_twosize[,null_cpns,drop = FALSE]) / rowSums(z_cpn_twosize)
  FDR_left = rowSums(z_cpn_ldist[,null_cpns,drop = FALSE]) / rowSums(z_cpn_ldist)
  FDR_right = rowSums(z_cpn_rdist[,null_cpns,drop = FALSE]) / rowSums(z_cpn_rdist)
  
  ## plot
  if (plot) {
    par(mfrow = c(1,3))
    zz = hist(z, breaks = 100)
    zz_cpn_prob = matrix(0, nrow = length(zz$mids), ncol = J)
    for (j in 1:J) zz_cpn_prob[,j] = piprob[j] * dnorm(zz$mids, mean = mu[j], sd = sqrt(sigma2[j] + 1))
    lfdr_zz = rowSums(zz_cpn_prob[,null_cpns,drop = FALSE]) / rowSums(zz_cpn_prob)
    # plot(zz)
    for (i in 1:length(zz$mids)) lines(c(zz$mids[i], zz$mids[i]), c(0, (1 - lfdr_zz[i]) * zz$counts[i]), col = 6)
    
    plot(sort(z), lfdr[order(z)], type = 'l', lwd = 1)
    lines(sort(z), FDR[order(z)], type = 'l', lty = 2, lwd = 1)
    lines(sort(z), FDR_left[order(z)], lty = 2, lwd = 1, col = "red")
    lines(sort(z), FDR_right[order(z)], lty = 2, lwd = 1, col = "blue")
    legend("topleft", legend = c("lfdr", "two side", "left", "right"), col = c("black", "black", "red", "blue"), 
          lty = c(1,2,2,2), lwd = 1, bty = "n")
    
    hist(effect_sz[order(z)], breaks = 100)
  }
  
  res = list(lfdr = lfdr, effect_sz = effect_sz, var_effect_sz = var_effect_sz, 
             FDR = FDR, FDR_left = FDR_left, FDR_right = FDR_right,
             prior_density = mix_prior_den, mix_prior = mix_prior, null_cpns = null_cpns)
  if (turning_P) res[["turning_P"]] = turn_res
  return(res)
}

