suppressWarnings(suppressPackageStartupMessages(library(Rcpp)))
suppressWarnings(suppressPackageStartupMessages(library(RcppArmadillo)))

source('utilities.R')
Rcpp::sourceCpp('utilities.cpp')

HOM_HOV <- function(
    dat,
    fixed_eff_ind,
    rand_eff_ind,
    n_sim,
    n_burn
) {
  eta_u <- eta_w <- eta_z <- 3
  w_perc <- 3
  is_outlier_only <- FALSE
  
  mod <- EVIII_RE_sampler_II(
    dat, NULL, 
    fixed_eff_ind, rand_eff_ind, 
    1,1,1,
    # fixed_eff_u_ind, fixed_eff_w_ind, fixed_eff_z_ind,
    n_sim, 10000, 123456, 
    diag(1,length(rand_eff_ind)),
    eta_u, eta_w, eta_z, 
    FALSE, NULL,
    n_sim, n_sim, n_sim, 
    w_perc, is_outlier_only
  )
  
  n <- length(dat$time)
  ni_arr <- sapply(dat$time, length)
  ni_max <- max(ni_arr)
  na_cols_list <- list()
  for (i in 1:n) {
    if (ni_arr[i]<ni_max) { na_cols_list[[i]] <- ((i-1)*ni_max) + (ni_arr[i]+1):ni_max }
  }
  na_cols <- do.call(c, na_cols_list)
  mod$mcmc.df$w[,na_cols] <- NA
  
  idx <- (n_burn+1):n_sim
  return(list(
    alpha = mod$mcmc.df$alpha[idx,],
    beta = mod$mcmc.df$beta[idx,],
    gamma.z = mod$mcmc.df$gamma.z[idx,,drop=F],
    gamma.u = mod$mcmc.df$gamma.u[idx,,drop=F],
    gamma.w = mod$mcmc.df$gamma.w[idx,,drop=F],
    Lambda = mod$mcmc.df$Lambda[idx,],
    r_mat = mod$mcmc.df$r_mat[idx,],
    rho = mod$mcmc.df$rho[idx],
    sigma.e2 = mod$mcmc.df$sigma.e2[idx,],
    sigma2 = mod$mcmc.df$sigma2[idx,],
    u = mod$mcmc.df$u[idx,],
    w = mod$mcmc.df$w[idx,],
    z = mod$mcmc.df$z[idx,]
  ))
}


HEM_HOV <- function(
    dat,
    fixed_eff_ind,
    rand_eff_ind,
    fixed_eff_u_ind,
    eta_u,
    n_sim,
    n_burn
) {
  eta_w <- eta_z <- 3 ## setting to some default values, but will not be used while sampling
  w_perc <- 3
  is_outlier_only <- FALSE
  
  mod <- EVIII_RE_sampler_II(
    dat, NULL, 
    fixed_eff_ind, rand_eff_ind, 
    fixed_eff_u_ind, 1, 1,
    n_sim, 10000, 123456, 
    diag(1,length(rand_eff_ind)),
    eta_u, eta_w, eta_z, 
    FALSE, NULL,
    0, n_sim, n_sim, 
    w_perc, is_outlier_only
  )
  
  n <- length(dat$time)
  ni_arr <- sapply(dat$time, length)
  ni_max <- max(ni_arr)
  na_cols_list <- list()
  for (i in 1:n) {
    if (ni_arr[i]<ni_max) { na_cols_list[[i]] <- ((i-1)*ni_max) + (ni_arr[i]+1):ni_max }
  }
  na_cols <- do.call(c, na_cols_list)
  mod$mcmc.df$w[,na_cols] <- NA
  
  idx <- (n_burn+1):n_sim
  return(list(
    alpha = mod$mcmc.df$alpha[idx,],
    beta = mod$mcmc.df$beta[idx,],
    gamma.z = mod$mcmc.df$gamma.z[idx,,drop=F],
    gamma.u = mod$mcmc.df$gamma.u[idx,,drop=F],
    gamma.w = mod$mcmc.df$gamma.w[idx,,drop=F],
    Lambda = mod$mcmc.df$Lambda[idx,],
    r_mat = mod$mcmc.df$r_mat[idx,],
    rho = mod$mcmc.df$rho[idx],
    sigma.e2 = mod$mcmc.df$sigma.e2[idx,],
    sigma2 = mod$mcmc.df$sigma2[idx,],
    u = mod$mcmc.df$u[idx,],
    w = mod$mcmc.df$w[idx,],
    z = mod$mcmc.df$z[idx,]
  ))
}


HOM_HOV_O <- function(
    dat,
    fixed_eff_ind,
    rand_eff_ind,
    fixed_eff_w_ind,
    eta_w,
    n_sim,
    n_burn
) {
  eta_u <- eta_z <- 3 ## setting to some default values, but will not be used while sampling
  w_perc <- 3
  is_outlier_only <- TRUE
  
  mod <- EVIII_RE_sampler_II(
    dat, NULL, 
    fixed_eff_ind, rand_eff_ind, 
    1, fixed_eff_w_ind, 1,
    n_sim, 10000, 123456, 
    diag(1,length(rand_eff_ind)),
    eta_u, eta_w, eta_z, 
    FALSE, NULL,
    0, n_sim, n_sim, 
    w_perc, is_outlier_only
  )
  
  n <- length(dat$time)
  ni_arr <- sapply(dat$time, length)
  ni_max <- max(ni_arr)
  na_cols_list <- list()
  for (i in 1:n) {
    if (ni_arr[i]<ni_max) { na_cols_list[[i]] <- ((i-1)*ni_max) + (ni_arr[i]+1):ni_max }
  }
  na_cols <- do.call(c, na_cols_list)
  mod$mcmc.df$w[,na_cols] <- NA
  
  idx <- (n_burn+1):n_sim
  return(list(
    alpha = mod$mcmc.df$alpha[idx,],
    beta = mod$mcmc.df$beta[idx,],
    gamma.z = mod$mcmc.df$gamma.z[idx,,drop=F],
    gamma.u = mod$mcmc.df$gamma.u[idx,,drop=F],
    gamma.w = mod$mcmc.df$gamma.w[idx,,drop=F],
    Lambda = mod$mcmc.df$Lambda[idx,],
    r_mat = mod$mcmc.df$r_mat[idx,],
    rho = mod$mcmc.df$rho[idx],
    sigma.e2 = mod$mcmc.df$sigma.e2[idx,],
    sigma2 = mod$mcmc.df$sigma2[idx,],
    u = mod$mcmc.df$u[idx,],
    w = mod$mcmc.df$w[idx,],
    z = mod$mcmc.df$z[idx,]
  ))
}


HEM_HOV_O <- function(
    dat,
    fixed_eff_ind,
    rand_eff_ind,
    fixed_eff_u_ind,
    fixed_eff_w_ind,
    eta_u,
    eta_w,
    n_sim,
    n_burn
) {
  eta_z <- 3 ## setting to some default values, but will not be used while sampling
  w_perc <- 3
  is_outlier_only <- TRUE
  
  mod <- EVIII_RE_sampler_II(
    dat, NULL, 
    fixed_eff_ind, rand_eff_ind, 
    fixed_eff_u_ind, fixed_eff_w_ind, 1,
    n_sim, 10000, 123456, 
    diag(1,length(rand_eff_ind)),
    eta_u, eta_w, eta_z, 
    FALSE, NULL,
    0, 0, n_sim,
    w_perc, is_outlier_only
  )
  
  n <- length(dat$time)
  ni_arr <- sapply(dat$time, length)
  ni_max <- max(ni_arr)
  na_cols_list <- list()
  for (i in 1:n) {
    if (ni_arr[i]<ni_max) { na_cols_list[[i]] <- ((i-1)*ni_max) + (ni_arr[i]+1):ni_max }
  }
  na_cols <- do.call(c, na_cols_list)
  mod$mcmc.df$w[,na_cols] <- NA
  
  idx <- (n_burn+1):n_sim
  return(list(
    alpha = mod$mcmc.df$alpha[idx,],
    beta = mod$mcmc.df$beta[idx,],
    gamma.z = mod$mcmc.df$gamma.z[idx,,drop=F],
    gamma.u = mod$mcmc.df$gamma.u[idx,,drop=F],
    gamma.w = mod$mcmc.df$gamma.w[idx,,drop=F],
    Lambda = mod$mcmc.df$Lambda[idx,],
    r_mat = mod$mcmc.df$r_mat[idx,],
    rho = mod$mcmc.df$rho[idx],
    sigma.e2 = mod$mcmc.df$sigma.e2[idx,],
    sigma2 = mod$mcmc.df$sigma2[idx,],
    u = mod$mcmc.df$u[idx,],
    w = mod$mcmc.df$w[idx,],
    z = mod$mcmc.df$z[idx,]
  ))
}


HEM_HEV_O <- function(
    dat,
    fixed_eff_ind,
    rand_eff_ind,
    fixed_eff_u_ind,
    fixed_eff_w_ind,
    fixed_eff_z_ind,
    eta_u,
    eta_w,
    eta_z,
    n_sim,
    n_burn
) {
  w_perc <- 3
  is_outlier_only <- TRUE
  
  mod <- EVIII_RE_sampler_II(
    dat, NULL, 
    fixed_eff_ind, rand_eff_ind, 
    fixed_eff_u_ind, fixed_eff_w_ind, fixed_eff_z_ind,
    n_sim, 10000, 123456, 
    diag(1,length(rand_eff_ind)),
    eta_u, eta_w, eta_z, 
    FALSE, NULL,
    0, 0, 0,
    w_perc, is_outlier_only
  )
  
  n <- length(dat$time)
  ni_arr <- sapply(dat$time, length)
  ni_max <- max(ni_arr)
  na_cols_list <- list()
  for (i in 1:n) {
    if (ni_arr[i]<ni_max) { na_cols_list[[i]] <- ((i-1)*ni_max) + (ni_arr[i]+1):ni_max }
  }
  na_cols <- do.call(c, na_cols_list)
  mod$mcmc.df$w[,na_cols] <- NA
  
  idx <- (n_burn+1):n_sim
  return(list(
    alpha = mod$mcmc.df$alpha[idx,],
    beta = mod$mcmc.df$beta[idx,],
    gamma.z = mod$mcmc.df$gamma.z[idx,,drop=F],
    gamma.u = mod$mcmc.df$gamma.u[idx,,drop=F],
    gamma.w = mod$mcmc.df$gamma.w[idx,,drop=F],
    Lambda = mod$mcmc.df$Lambda[idx,],
    r_mat = mod$mcmc.df$r_mat[idx,],
    rho = mod$mcmc.df$rho[idx],
    sigma.e2 = mod$mcmc.df$sigma.e2[idx,],
    sigma2 = mod$mcmc.df$sigma2[idx,],
    u = mod$mcmc.df$u[idx,],
    w = mod$mcmc.df$w[idx,],
    z = mod$mcmc.df$z[idx,]
  ))
}


extract_params <- function (mod, mod_type) {
  if (mod_type == "HOM-HOV") {
    params_list <- list(
      alpha = mod$alpha[,2],
      beta = mod$beta,
      Lambda = mod$Lambda,
      r_mat = mod$r_mat,
      rho = mod$beta,
      sigma02 = mod$sigma2[,1]
    )
  } else if (mod_type == "HEM-HOV") {
    params_list <- list(
      alpha = mod$alpha[,2],
      beta = mod$beta,
      gamma.u = mod$gamma.u,
      Lambda = mod$mcmc.df$Lambda,
      r_mat = mod$mcmc.df$r_mat,
      rho = mod$mcmc.df$beta,
      sigma02 = mod$mcmc.df$sigma2[,1],
      u = mod$mcmc.df$u
    )
  } else if (mod_type == "HOM-HOV-O") {
    params_list <- list(
      alpha = mod$mcmc.df$alpha[,2],
      beta = mod$mcmc.df$beta,
      gamma.w = mod$mcmc.df$gamma.w,
      Lambda = mod$mcmc.df$Lambda,
      r_mat = mod$mcmc.df$r_mat,
      rho = mod$mcmc.df$beta,
      sigma02 = mod$mcmc.df$sigma2[,1],
      w = mod$mcmc.df$w
    )
  } else if (mod_type == "HEM-HOV-O") {
    params_list <- list(
      alpha = mod$mcmc.df$alpha[,2],
      beta = mod$mcmc.df$beta,
      gamma.u = mod$mcmc.df$gamma.u,
      gamma.w = mod$mcmc.df$gamma.w,
      Lambda = mod$mcmc.df$Lambda,
      r_mat = mod$mcmc.df$r_mat,
      rho = mod$mcmc.df$beta,
      sigma02 = mod$mcmc.df$sigma2[,1,drop=F],
      u = mod$mcmc.df$u,
      w = mod$mcmc.df$w
    )
  } else if (mod_type == "HEM-HEV-O") {
    params_list <- list(
      alpha = mod$mcmc.df$alpha[,2],
      beta = mod$mcmc.df$beta,
      gamma.u = mod$mcmc.df$gamma.u,
      gamma.w = mod$mcmc.df$gamma.w,
      gamma.z = mod$mcmc.df$gamma.h,
      Lambda = mod$mcmc.df$Lambda,
      r_mat = mod$mcmc.df$r_mat,
      rho = mod$mcmc.df$beta,
      sigma02 = mod$mcmc.df$sigma2[,1],
      sigma12 = mod$mcmc.df$sigma2[,2],
      u = mod$mcmc.df$u,
      w = mod$mcmc.df$w,
      z = mod$mcmc.df$z
    )
  }
}


compute_WAIC_HOM_HOV <- function(
    dat,
    mod,
    iter_arr,
    fixed_eff_ind,
    rand_eff_ind,
    type = "C-WAIC"
) {
  
  n <- length(dat$data)
  ni <- sapply(dat$data, length)
  ni_max <- max(ni)
  eta_u <- eta_w <- eta_z <- 3
  mod_type <- 1
  if (type == "C-WAIC") {
    is_conditional <- 1
  } else if (type == "M-WAIC") {
    is_conditional <- 0
  }
  
  W_list <- list()
  for (i in seq_along(ni)) {
    mat <- matrix(0, nrow = ni[i]*(ni[i]+1)/2, ncol = ni[i])
    count <- 1
    for (j in 1:ni[i]) {
      for (k in 1:j) {
        mat[count,c(j,k)] <- 1
        count <- count + 1
      }
    }
    W_list[[i]] <- c(
      0,
      sapply(
        1:nrow(mat),
        function(ind) {
          sum(mat[ind,] * 2^((ni[i]-1):0))
        }
      )
    )
  }
  
  WAIC <- compute_WAIC_EVIII_RE_marginal_1(
    dat, mod, W_list, iter_arr,
    fixed_eff_ind, rand_eff_ind,
    1, 1, 1,
    eta_w, eta_u, eta_z, mod_type, is_conditional
  )
  return(WAIC)
}



compute_WAIC_HEM_HOV <- function(
    dat,
    mod,
    iter_arr,
    fixed_eff_ind,
    rand_eff_ind,
    fixed_eff_u_ind,
    eta_u,
    type = "C-WAIC"
) {
  
  n <- length(dat$data)
  ni <- sapply(dat$data, length)
  ni_max <- max(ni)
  eta_w <- eta_z <- 3
  mod_type <- 2
  if (type == "C-WAIC") {
    is_conditional <- 1
  } else if (type == "M-WAIC") {
    is_conditional <- 0
  }
  
  W_list <- list()
  for (i in seq_along(ni)) {
    mat <- matrix(0, nrow = ni[i]*(ni[i]+1)/2, ncol = ni[i])
    count <- 1
    for (j in 1:ni[i]) {
      for (k in 1:j) {
        mat[count,c(j,k)] <- 1
        count <- count + 1
      }
    }
    W_list[[i]] <- c(
      0,
      sapply(
        1:nrow(mat),
        function(ind) {
          sum(mat[ind,] * 2^((ni[i]-1):0))
        }
      )
    )
  }
  
  WAIC <- compute_WAIC_EVIII_RE_marginal_1(
    dat, mod, W_list, iter_arr,
    fixed_eff_ind, rand_eff_ind,
    fixed_eff_u_ind, 1, 1,
    eta_w, eta_u, eta_z, mod_type, is_conditional
  )
  return(WAIC)
}


compute_WAIC_HOM_HOV_O <- function(
    dat,
    mod,
    iter_arr,
    fixed_eff_ind,
    rand_eff_ind,
    fixed_eff_w_ind,
    eta_w,
    type = "C-WAIC"
) {
  
  n <- length(dat$data)
  ni <- sapply(dat$data, length)
  ni_max <- max(ni)
  eta_u <- eta_z <- 3
  mod_type <- 5
  if (type == "C-WAIC") {
    is_conditional <- 1
  } else if (type == "M-WAIC") {
    is_conditional <- 0
  }
  
  W_list <- list()
  for (i in seq_along(ni)) {
    mat <- matrix(0, nrow = ni[i]*(ni[i]+1)/2, ncol = ni[i])
    count <- 1
    for (j in 1:ni[i]) {
      for (k in 1:j) {
        mat[count,c(j,k)] <- 1
        count <- count + 1
      }
    }
    W_list[[i]] <- c(
      0,
      sapply(
        1:nrow(mat),
        function(ind) {
          sum(mat[ind,] * 2^((ni[i]-1):0))
        }
      )
    )
  }
  
  WAIC <- compute_WAIC_EVIII_RE_marginal_1(
    dat, mod, W_list, iter_arr,
    fixed_eff_ind, rand_eff_ind,
    1, fixed_eff_w_ind, 1,
    eta_w, eta_u, eta_z, mod_type, is_conditional
  )
  return(WAIC)
}


compute_WAIC_HEM_HOV_O <- function(
    dat,
    mod,
    iter_arr,
    fixed_eff_ind,
    rand_eff_ind,
    fixed_eff_u_ind,
    fixed_eff_w_ind,
    eta_u,
    eta_w,
    type = "C-WAIC"
) {
  
  n <- length(dat$data)
  ni <- sapply(dat$data, length)
  ni_max <- max(ni)
  eta_z <- 3
  mod_type <- 3
  if (type == "C-WAIC") {
    is_conditional <- 1
  } else if (type == "M-WAIC") {
    is_conditional <- 0
  }
  
  W_list <- list()
  for (i in seq_along(ni)) {
    mat <- matrix(0, nrow = ni[i]*(ni[i]+1)/2, ncol = ni[i])
    count <- 1
    for (j in 1:ni[i]) {
      for (k in 1:j) {
        mat[count,c(j,k)] <- 1
        count <- count + 1
      }
    }
    W_list[[i]] <- c(
      0,
      sapply(
        1:nrow(mat),
        function(ind) {
          sum(mat[ind,] * 2^((ni[i]-1):0))
        }
      )
    )
  }
  
  WAIC <- compute_WAIC_EVIII_RE_marginal_1(
    dat, mod, W_list, iter_arr,
    fixed_eff_ind, rand_eff_ind,
    fixed_eff_u_ind, fixed_eff_w_ind, 1,
    eta_w, eta_u, eta_z, mod_type, is_conditional
  )
  return(WAIC)
}


compute_WAIC_HEM_HEV_O <- function(
    dat,
    mod,
    iter_arr,
    fixed_eff_ind,
    rand_eff_ind,
    fixed_eff_u_ind,
    fixed_eff_w_ind,
    fixed_eff_z_ind,
    eta_u,
    eta_w,
    eta_z,
    type = "C-WAIC"
) {
  
  n <- length(dat$data)
  ni <- sapply(dat$data, length)
  ni_max <- max(ni)
  mod_type <- 4
  if (type == "C-WAIC") {
    is_conditional <- 1
  } else if (type == "M-WAIC") {
    is_conditional <- 0
  }
  
  W_list <- list()
  for (i in seq_along(ni)) {
    mat <- matrix(0, nrow = ni[i]*(ni[i]+1)/2, ncol = ni[i])
    count <- 1
    for (j in 1:ni[i]) {
      for (k in 1:j) {
        mat[count,c(j,k)] <- 1
        count <- count + 1
      }
    }
    W_list[[i]] <- c(
      0,
      sapply(
        1:nrow(mat),
        function(ind) {
          sum(mat[ind,] * 2^((ni[i]-1):0))
        }
      )
    )
  }
  
  WAIC <- compute_WAIC_EVIII_RE_marginal_1(
    dat, mod, W_list, iter_arr,
    fixed_eff_ind, rand_eff_ind,
    fixed_eff_u_ind, fixed_eff_w_ind, fixed_eff_z_ind,
    eta_w, eta_u, eta_z, mod_type, is_conditional
  )
  return(WAIC)
}

