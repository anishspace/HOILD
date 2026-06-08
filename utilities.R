library(MASS)
suppressWarnings(library(mvtnorm))
suppressWarnings(library(invgamma))
suppressWarnings(library(ggplot2))
suppressWarnings(library(ggpubr))
suppressWarnings(library(gridExtra))
suppressWarnings(library(ggExtra))

view.doc <- function (file.path, type='html') {
  file.dir.path <- dirname(file.path)
  file.name <- basename(file.path)
  options(knitr.duplicate.label = 'allow')
  rmarkdown::render(file.path, output_format = paste(type,'document',sep='_'))
  doc.file.name <- paste(strsplit(file.name, '.', fixed=TRUE)[[1]][1],
                         type, sep='.')
  doc.file.path <- paste(file.dir.path, doc.file.name, sep='/')
  if (file.copy(doc.file.path, file.path(tempdir(), doc.file.name), overwrite=TRUE)) {
    viewer <- getOption('viewer')
    viewer(file.path(tempdir(), doc.file.name))
  } else {
    stop ("error: failed to copy file to the tmp directory.\n")
  }
  
}

get.rho.est <- function(data) {
  rho <- array(0, dim(data)[1])
  for (i in 1:dim(data)[1]) {
    rho[i] <- acf(data[i,], plot=FALSE)$acf[2]
  }
  return (rho)
}

H.mat <- function(ni, tt, sd, rho) {
  Hmat <- matrix(0, nrow=ni, ncol=ni)
  for (i in 1:ni) {
    for (j in 1:ni) {
      Hmat[i,j] = sd^2 * rho^(abs(tt[j]-tt[i]))
    }
  }
  return (Hmat)
}

## drop.rate will be #missing/#total
## to have 50% dropout, we should ni.min from 1 to 10, so that on average
## there are 5 observations per person, i.e., 50% missing
sim.data <- function(n, ni.max, drop.rate, vary.ni, beta, asc1, asc2, z,
                     sd.serial.vec, is.exact=TRUE, rho.serial.vec, prob.serial = NULL, 
                     sd.error, plot=TRUE) {
  if (vary.ni) {
    ## When number of missing observations is also random, with overall drop rate
    ## equal to drop.rate
    ni.min <- 1 + max(0, ceiling((ni.max-1)*(1-2*drop.rate)))
    if (ni.min == ni.max) {
      ni <- sample(c(ni.max,ni.max), n, replace=TRUE)
    } else {
      ni <- sample(seq(ni.min, ni.max, 1), n, replace=TRUE)
    }
  } else {
    # ni is fixed for all individuals, that is #missing observations is fixed
    ni <- rep(max(1, ceiling(ni.max*(1-drop.rate))), n)
  }
  
  # when number of missing observation is fixed and equal to ni.max*droprate
  tt <- lapply(1:n, function(i) sort(sample(1:ni.max,ni[i],replace=FALSE)))
  if (is.exact) {
    serial.sd <- sd.serial.vec
  } else {
    serial.sd <- sample(sd.serial.vec, n, replace=TRUE, prob.serial)
  }
  serial.rho <- sample(rho.serial.vec, n, replace=TRUE) ## Unform distribution of rho
  serial.eff <- matrix(NA, nrow=n, ncol=ni.max)
  
  # fixed effects
  fixed.eff <- list()
  for (i in 1:n) {
    fixed.eff[[i]] <- matrix(NA, nrow=ni.max, ncol=5)
    fixed.eff[[i]][tt[[i]],1] <- rep(1, ni[i])
    fixed.eff[[i]][tt[[i]],2] <- tt[[i]]
    fixed.eff[[i]][tt[[i]],3] <- tt[[i]]^2
    fixed.eff[[i]][tt[[i]],4] <- rep(sample(c(z[i],1-z[i]), 1, 
                                            prob=c(asc1, 1-asc1)), ni[i])
    fixed.eff[[i]][tt[[i]],5] <- rep(sample(c(z[i],1-z[i]), 1, 
                                            prob=c(asc2, 1-asc2)), ni[i])
  }
  
  dat <- matrix(NA, nrow=n, ncol=ni.max)
  err.sd <- sample(sd.error, n, replace=TRUE)
  for (i in 1:n) {
    serial.eff[i,tt[[i]]] <- mvrnorm(1, mu=rep(0,ni[i]), 
                                     Sigma=H.mat(ni[i], tt[[i]], 
                                                 serial.sd[i], serial.rho[i]))
    dat[i,tt[[i]]] <- mvrnorm(1, mu = fixed.eff[[i]][tt[[i]],]%*%beta + 
                                serial.eff[i,tt[[i]]], 
                              Sigma=diag(err.sd[i]^2, ni[i]))
  }
  
  dat.trend <- list(data = dat, time = tt, beta=beta, z = z, asc1 = asc1, asc2 = asc2,
                    fixed.eff = fixed.eff,  serial.eff = serial.eff, 
                    serial.eff.sd = serial.sd, serial.rho = serial.rho, err.sd = err.sd)
  # if (plot) diagnostic.data(dat.trend)
  
  return (dat.trend)
}

sim.error.data <- function(n, ni.max, drop.rate, vary.ni, beta, asc1, asc2, z,
                           sd.error.vec, rho.error.vec, prob.error = NULL) {
  if (vary.ni) {
    ## When number of missing observations is also random, with overall drop rate
    ## equal to drop.rate
    ni.min <- 1 + max(0, ceiling((ni.max-1)*(1-2*drop.rate)))
    if (ni.min == ni.max) {
      ni <- sample(c(ni.max,ni.max), n, replace=TRUE)
    } else {
      ni <- sample(seq(ni.min, ni.max, 1), n, replace=TRUE)
    }
  } else {
    # ni is fixed for all individuals, that is #missing observations is fixed
    ni <- rep(max(1, ceiling(ni.max*(1-drop.rate))), n)
  }
  
  # when number of missing observation is fixed and equal to ni.max*droprate
  tt <- lapply(1:n, function(i) sort(sample(1:ni.max,ni[i],replace=FALSE)))
  
  # fixed effects
  fixed.eff <- list()
  for (i in 1:n) {
    fixed.eff[[i]] <- matrix(NA, nrow=ni.max, ncol=5)
    fixed.eff[[i]][tt[[i]],1] <- rep(1, ni[i])
    fixed.eff[[i]][tt[[i]],2] <- tt[[i]]
    fixed.eff[[i]][tt[[i]],3] <- tt[[i]]^2
    fixed.eff[[i]][tt[[i]],4] <- rep(sample(c(z[i],1-z[i]), 1, 
                                            prob=c(asc1, 1-asc1)), ni[i])
    fixed.eff[[i]][tt[[i]],5] <- rep(sample(c(z[i],1-z[i]), 1, 
                                            prob=c(asc2, 1-asc2)), ni[i])
  }
  
  error.sd <- sample(sd.error.vec, n, replace=TRUE, prob.error)
  error.rho <- sample(rho.error.vec, n, replace=TRUE) ## Unform distribution of rho
  
  dat <- matrix(NA, nrow=n, ncol=ni.max)
  for (i in 1:n) {
    dat[i,tt[[i]]] <- mvrnorm(1, mu = fixed.eff[[i]][tt[[i]],]%*%beta, 
                              Sigma = H.mat(ni[i], tt[[i]], error.sd[i], error.rho[i]))
  }
  return (list(data = dat, time = tt, beta=beta, z = z, asc1 = asc1, asc2 = asc2, 
               fixed.eff = fixed.eff, error.sd = error.sd, error.rho = error.rho))
}

sim.error.dataIV <- function(n, ni.max, drop.rate, vary.ni, beta, fixed.eff.z.ind, gam,
                             sd.error, is.exact=TRUE, rho.error.vec, prob.error = NULL) {
  if (vary.ni) {
    ## When number of missing observations is also random, with overall drop rate
    ## equal to drop.rate
    ni.min <- 1 + max(0, ceiling((ni.max-1)*(1-2*drop.rate)))
    if (ni.min == ni.max) {
      ni <- sample(c(ni.max,ni.max), n, replace=TRUE)
    } else {
      ni <- sample(seq(ni.min, ni.max, 1), n, replace=TRUE)
    }
  } else {
    # ni is fixed for all individuals, that is #missing observations is fixed
    ni <- rep(max(1, ceiling(ni.max*(1-drop.rate))), n)
  }
  
  # when number of missing observation is fixed and equal to ni.max*droprate
  tt <- lapply(1:n, function(i) sort(sample(1:ni.max,ni[i],replace=FALSE)))
  
  # fixed effects
  fixed.eff <- list()
  z <- array(NA,n)
  sd.error.vec <- array(NA,n)
  for (i in 1:n) {
    fixed.eff[[i]] <- matrix(NA, nrow=ni.max, ncol=5)
    fixed.eff[[i]][tt[[i]],1] <- rep(1, ni[i])
    fixed.eff[[i]][tt[[i]],2] <- tt[[i]]
    fixed.eff[[i]][tt[[i]],3] <- tt[[i]]^2
    fixed.eff[[i]][tt[[i]],4] <- rep(sample(c(0,1), 1), ni[i])
    fixed.eff[[i]][tt[[i]],5] <- rep(sample(c(0,1), 1), ni[i])
    pi.p <- exp(fixed.eff[[i]][1,fixed.eff.z.ind] %*% gam) /
      (1 + exp(fixed.eff[[i]][1,fixed.eff.z.ind] %*% gam))
    z[i] <- sample(c(1,0), 1, prob=c(pi.p, 1-pi.p))
    sd.error.vec[i] <- ifelse(z[i]==1, 
                              sd.error[1], 
                              sqrt(rgamma(1, 1, 1/(sd.error[2]^2))))
  }
  
  if (is.exact) {
    error.sd <- sd.error.vec
  } else {
    error.sd <- sample(sd.error.vec, n, replace=TRUE, prob.error)
  }
  error.rho <- sample(rho.error.vec, n, replace=TRUE) ## Unform distribution of rho
  
  dat <- matrix(NA, nrow=n, ncol=ni.max)
  for (i in 1:n) {
    dat[i,tt[[i]]] <- mvrnorm(1, mu = fixed.eff[[i]][tt[[i]],]%*%beta, 
                              Sigma = H.mat(ni[i], tt[[i]], error.sd[i], error.rho[i]))
  }
  return (list(data = dat, time = tt, beta=beta, z = z, fixed.eff = fixed.eff, 
               error.sd = error.sd, error.rho = error.rho))
}

sim.error.dataV <- function(n, ni.max, drop.rate, vary.ni, pi.e, beta, fixed.eff.z.ind, 
                            gam, sd.error, is.exact=TRUE, rho.error.vec, 
                            prob.error = NULL) {
  set.seed(123456)
  if (vary.ni) {
    ## When number of missing observations is also random, with overall drop rate
    ## equal to drop.rate
    ni.min <- 1 + max(0, ceiling((ni.max-1)*(1-2*drop.rate)))
    if (ni.min == ni.max) {
      ni <- sample(c(ni.max,ni.max), n, replace=TRUE)
    } else {
      ni <- sample(seq(ni.min, ni.max, 1), n, replace=TRUE)
    }
  } else {
    # ni is fixed for all individuals, that is #missing observations is fixed
    ni <- rep(max(1, ceiling(ni.max*(1-drop.rate))), n)
  }
  
  # when number of missing observation is fixed and equal to ni.max*droprate
  tt <- lapply(1:n, function(i) sort(sample(1:ni.max,ni[i],replace=FALSE)))
  
  # fixed effects
  fixed.eff <- list()
  sd.error.vec <- array(NA, n)
  z <- array(NA,n)
  pi.h <- array(NA,n)
  for (i in 1:n) {
    fixed.eff[[i]] <- matrix(NA, nrow=ni[i], ncol=5)
    fixed.eff[[i]][,1] <- rep(1, ni[i])
    fixed.eff[[i]][,2] <- tt[[i]]
    fixed.eff[[i]][,3] <- tt[[i]]^2
    fixed.eff[[i]][,4] <- rep(sample(c(0,1), 1), ni[i])
    fixed.eff[[i]][,5] <- rep(sample(c(0,1), 1), ni[i])
    pi.h <- exp(fixed.eff[[i]][1,fixed.eff.z.ind] %*% gam) /
      (1 + exp(fixed.eff[[i]][1,fixed.eff.z.ind] %*% gam))
    z[i] <- sample(c(1,0), 1, prob=c(pi.h, 1-pi.h))
    sd.error.vec[i] <- ifelse(z[i]==1, 
                              sd.error[1], 
                              sqrt(rgamma(1, 1, 1/(sd.error[2]^2))))
  }
  
  if (is.exact) {
    error.sd <- sd.error.vec
  } else {
    error.sd <- sample(sd.error.vec, n, replace=TRUE, prob.error)
  }
  error.rho <- sample(rho.error.vec, n, replace=TRUE) ## Unform distribution of rho
  
  # dat <- matrix(NA, nrow=n, ncol=ni.max)
  dat <- list()
  # cov.type <- sample(c(0,1), n, replace=TRUE, prob=c(1-pi.e, pi.e))
  # for (i in 1:n) {
  #   if (cov.type[i]) {
  #     dat[i,tt[[i]]] <- mvrnorm(1, mu = fixed.eff[[i]][tt[[i]],]%*%beta,
  #                               Sigma = diag(error.sd[i]^2, ni[i]))
  #   } else {
  #     dat[i,tt[[i]]] <- mvrnorm(1, mu = fixed.eff[[i]][tt[[i]],]%*%beta,
  #                               Sigma = H.mat(ni[i], tt[[i]], 
  #                                             error.sd[i], error.rho[i]))
  #   }
  # }
  for (i in 1:n) {
    # dat[i,tt[[i]]] <- mvrnorm(1, mu = fixed.eff[[i]][tt[[i]],]%*%beta, 
    #                           Sigma = pi.e*diag(error.sd[i]^2, ni[i]) +
    #                             (1-pi.e)*H.mat(ni[i], tt[[i]], error.sd[i], error.rho[i]))
    dat[[i]] <- mvrnorm(1, mu = fixed.eff[[i]]%*%beta, 
                        Sigma = pi.e*diag(error.sd[i]^2, ni[i]) +
                          (1-pi.e)*H.mat(ni[i], tt[[i]], error.sd[i], error.rho[i]))
  }
  return (list(data = dat, time = tt, beta=beta, pi.e=pi.e, z = z, fixed.eff = fixed.eff, 
               gam = gam, error.sd = error.sd, error.rho = error.rho))
}

sim.mixed.data <- function(n, ni.max, drop.rate, vary.ni, beta, 
                           sd.error, is.exact=TRUE, rho.error.vec, 
                           prob.error = NULL) {
  set.seed(123456)
  if (vary.ni) {
    ## When number of missing observations is also random, with overall drop rate
    ## equal to drop.rate
    ni.min <- 1 + max(0, ceiling((ni.max-1)*(1-2*drop.rate)))
    if (ni.min == ni.max) {
      ni <- sample(c(ni.max,ni.max), n, replace=TRUE)
    } else {
      ni <- sample(seq(ni.min, ni.max, 1), n, replace=TRUE)
    }
  } else {
    # ni is fixed for all individuals, that is #missing observations is fixed
    ni <- rep(max(1, ceiling(ni.max*(1-drop.rate))), n)
  }
  
  # when number of missing observation is fixed and equal to ni.max*droprate
  tt <- lapply(1:n, function(i) sort(sample(1:ni.max,ni[i],replace=FALSE)))
  
  # fixed effects
  fixed.eff <- list()
  sd.error.vec <- array(NA, n)
  z <- array(NA,n)
  pi.h <- array(NA,n)
  for (i in 1:n) {
    fixed.eff[[i]] <- matrix(NA, nrow=ni[i], ncol=5)
    fixed.eff[[i]][,1] <- rep(1, ni[i])
    fixed.eff[[i]][,2] <- tt[[i]]
    fixed.eff[[i]][,3] <- tt[[i]]^2
    fixed.eff[[i]][,4] <- rep(sample(c(0,1), 1), ni[i])
    fixed.eff[[i]][,5] <- rep(sample(c(0,1), 1), ni[i])
  }
  z <- sample(1:length(sd.error), n, replace=TRUE, prob=NULL)
  sd.error.vec <- sd.error[z]
  
  if (is.exact) {
    error.sd <- sd.error.vec
  } else {
    error.sd <- sample(sd.error.vec, n, replace=TRUE, prob.error)
  }
  error.rho <- sample(rho.error.vec, n, replace=TRUE) ## Unform distribution of rho
  
  # dat <- matrix(NA, nrow=n, ncol=ni.max)
  dat <- list()
  for (i in 1:n) {
    # dat[i,tt[[i]]] <- mvrnorm(1, mu = fixed.eff[[i]][tt[[i]],]%*%beta, 
    #                           Sigma = H.mat(ni[i], tt[[i]], error.sd[i], error.rho[i]))
    dat[[i]] <- mvrnorm(1, mu = fixed.eff[[i]]%*%beta, 
                        Sigma = H.mat(ni[i], tt[[i]], error.sd[i], error.rho[i]))
  }
  return (list(data = dat, time = tt, beta=beta, 
               z = z, fixed.eff = fixed.eff, 
               error.sd = error.sd, error.rho = error.rho))
  # noise <- mvrnorm(n, rep(0, ni.max), Sigma = diag(100, ni.max))
  # wt <- matrix(sample(c(0,1), size = n*ni.max, replace=TRUE, prob=c(0.9, 0.1)),
  #              nrow = n, ncol=ni.max)
  # noise <- noise * wt
  # return (list(data = dat + noise, time = tt, beta=beta, 
  #              z = z, fixed.eff = fixed.eff, 
  #              error.sd = error.sd, error.rho = error.rho))
}

read_params <- function (
    fname, type = 1
) {
  df <- data.matrix(read.csv(fname)[,type+1])
  true_params <- list()
  true_params$beta <- df[1:5]
  true_params$alpha <- df[6:7]
  true_params$sigma_02 <- df[8]
  true_params$rho <- df[9]
  true_params$Lambda <- matrix(
    c(df[10], df[11], df[11], df[12]), nrow = 2
  )
  true_params$gamma_u <- df[13:15]
  true_params$gamma_w <- df[16:19]
  true_params$gamma_z <- df[20:22]
  true_params$eta_u <- df[23]
  true_params$eta_w <- df[24]
  true_params$eta_z <- df[25]
  return(true_params)
}

sim_data <- function(
    n, ni, p_miss, true_params, 
    fixed_eff_ind, rand_eff_ind, fixed_eff_u_ind, fixed_eff_w_ind, fixed_eff_z_ind,
    delta
) {
  dat <- list()
  dat$time <- list()
  dat$fixed.eff <- list()
  dat$data <- list()
  dat$r <- list()
  dat$eps <- list()
  dat$sigma_i2 <- array(NA, n)
  dat$u <- array(NA, n)
  dat$w <- list()
  dat$z <- array(NA, n)
  dat$pi_w <- list()
  pres_id_list <- list()
  for (i in 1:n) {
    miss_id <- sample(c(TRUE,FALSE), ni, replace=TRUE, prob=c(p_miss, 1-p_miss))
    pres_id_list[[i]] <- (1:ni)[!miss_id]
    dat$fixed.eff[[i]] <- cbind(
      1,
      sample(c(0,1),1),
      rnorm(1),
      1:ni,
      (1:ni)^2
    )
  }
  dat$fixed.eff[[i]] <- dat$fixed.eff[[i]][pres_id_list[[i]],]
  fixed_eff_mat <- do.call(rbind, dat$fixed.eff)
  fixed_eff_mat_sc <- cbind(
    1, 
    sapply(
      2:ncol(fixed_eff_mat),
      function(k) { scale(fixed_eff_mat[,k], center=TRUE, scale=TRUE) }
    )
  )
  pb <- txtProgressBar(1,n,style=3)
  for (i in 1:n) {
    dat$time[[i]] <- (1:ni)[pres_id_list[[i]]]
    dat$fixed.eff[[i]] <- fixed_eff_mat_sc[ni*(i-1)+pres_id_list[[i]], fixed_eff_ind]
    
    pi_u_i <- exp(dat$fixed.eff[[i]][1,fixed_eff_z_ind] %*% true_params$gamma_u)[1]
    pi_u_i <- pi_u_i/(1+pi_u_i)
    dat$u[i] <- sample(c(1,0), 1, prob=c(pi_u_i, 1-pi_u_i))
    # dat$u[i] <- 0
    pi_w_i <- exp(dat$fixed.eff[[i]][,fixed_eff_w_ind] %*% true_params$gamma_w)
    pi_w_i <- pi_w_i/(1+pi_w_i)
    dat$pi_w[[i]] <- pi_w_i
    flag <- TRUE
    while (flag) {
      dat$w[[i]] <- sapply(
        1:length(pres_id_list[[i]]),
        function(k) { sample(c(1,0), 1, prob=c(pi_w_i[k], 1-pi_w_i[k])) }
      )
      if (sum(dat$w[[i]]) <= 2) { flag <- FALSE }
    }
    # dat$w[[i]] <- rep(0, length(pres_id_list[[i]]))
    W_i_mat <- diag(dat$w[[i]])
    I_i <- diag(rep(1, length(pres_id_list[[i]])))
    I_W_i_mat <- I_i - W_i_mat
    pi_z_i <- exp(dat$fixed.eff[[i]][1,fixed_eff_z_ind] %*% true_params$gamma_z)[1]
    pi_z_i <- pi_z_i/(1+pi_z_i)
    dat$z[i] <- sample(c(1,0), 1, prob=c(pi_z_i, 1-pi_z_i))
    # dat$z[i] <- 0
    
    ## sampling sigma_i2
    if (dat$z[i] == 0) {
      dat$sigma_i2[i] <- true_params$sigma_02
      # sigma_i2 <- runif(1, true_params$sigma_02-delta, true_params$sigma_02+delta)
    } else {
      dat$sigma_i2[i] <- runif(
        1,
        true_params$sigma_02*true_params$eta_z^2 - delta*true_params$eta_z^2,
        true_params$sigma_02*true_params$eta_z^2 + delta*true_params$eta_z^2
      )
    }
    
    mean_i <- dat$fixed.eff[[i]][,fixed_eff_ind] %*% true_params$beta
    rand_eff_i <- dat$fixed.eff[[i]][,rand_eff_ind]
    
    flag <- TRUE
    while (flag) {
      r_i <- rmvn_cpp(
        1, 
        rep(0, length(rand_eff_ind)), 
        # true_params$sigma_02 *
        dat$sigma_i2[i] *
          true_params$eta_u^(2*dat$u[i]) * true_params$Lambda,
        1, FALSE
      )[1,]
      if (dat$u[i] == 0) { flag <- FALSE; next }
      if (dat$u[i] == 1) {
        quad <- t(r_i) %*% 
          (solve(true_params$sigma_02 * true_params$Lambda)) %*%
          # (solve(dat$sigma_i2[i] * true_params$Lambda)) %*% 
          r_i
        if (quad > qgamma(0.999, shape=1, rate=1/2)) { flag <- FALSE }
      }
    }
    dat$r[[i]] <- r_i
    
    flag <- TRUE
    flag1 <- flag2 <- TRUE
    while (flag) {
      eps_i <- rmvn_cpp(
        1, 
        rep(0, length(pres_id_list[[i]])), 
        dat$sigma_i2[i] *
          ((true_params$eta_w-1)*W_i_mat + I_i) %*% 
          h_mat(dat$time[[i]], true_params$rho) %*%
          ((true_params$eta_w-1)*W_i_mat + I_i),
        1, FALSE
      )[1,]
      w_ind <- as.logical(dat$w[[i]])
      if (dat$z[i] == 0) { 
        # if (all(w_ind) == FALSE) { flag <- FALSE; next }
        if (all(abs(eps_i[w_ind]) >  qnorm(
          0.9975, 0, sqrt(true_params$sigma_02)
        ))) { flag <- FALSE }
      } else if (dat$z[i] == 1) {
        quad <- 
          t(eps_i)[!w_ind] %*% (solve(
            true_params$sigma_02 *
            ((true_params$eta_w-1)*W_i_mat + I_i)[!w_ind, !w_ind] %*% 
              h_mat(dat$time[[i]], true_params$rho)[!w_ind, !w_ind] %*%
              ((true_params$eta_w-1)*W_i_mat + I_i)[!w_ind, !w_ind]
          )) %*% 
          eps_i[!w_ind]
        if (quad > qgamma(0.995, shape=sum(!w_ind)/2, rate=1/2)) { flag1 <- FALSE }
        if (all(abs(eps_i[w_ind]) >  qnorm(
          0.9975, 0, sqrt(dat$sigma_i2[i])
        ))) { flag2 <- FALSE }
        
        if(!flag1 & !flag2) { flag <- FALSE }
      }
    }
    dat$eps[[i]] <- eps_i
    
    dat$data[[i]] <- mean_i + rand_eff_i %*% dat$r[[i]] + dat$eps[[i]]
    
    setTxtProgressBar(pb, i)
  }
  dat$r_mat <- do.call(rbind, dat$r)
  return (dat)
}


sim_data1 <- function(
    n, ni, p_miss, true_params, 
    fixed_eff_ind, rand_eff_ind, fixed_eff_u_ind, fixed_eff_w_ind, fixed_eff_z_ind,
    delta
) {
  dat <- list()
  dat$time <- list()
  dat$fixed.eff <- list()
  dat$data <- list()
  dat$r <- list()
  dat$eps <- list()
  dat$sigma_i2 <- array(NA, n)
  dat$u <- array(NA, n)
  dat$w <- list()
  dat$z <- array(NA, n)
  dat$pi_w <- list()
  pres_id_list <- list()
  for (i in 1:n) {
    miss_id <- sample(c(TRUE,FALSE), ni, replace=TRUE, prob=c(p_miss, 1-p_miss))
    pres_id_list[[i]] <- (1:ni)[!miss_id]
    dat$fixed.eff[[i]] <- cbind(
      1,
      sample(c(0,1),1),
      rnorm(1),
      1:ni,
      (1:ni)^2
    )
  }
  dat$fixed.eff[[i]] <- dat$fixed.eff[[i]][pres_id_list[[i]],]
  fixed_eff_mat <- do.call(rbind, dat$fixed.eff)
  fixed_eff_mat_sc <- cbind(
    1, 
    sapply(
      2:ncol(fixed_eff_mat),
      function(k) { scale(fixed_eff_mat[,k], center=TRUE, scale=TRUE) }
    )
  )
  pb <- txtProgressBar(1,n,style=3)
  for (i in 1:n) {
    dat$time[[i]] <- (1:ni)[pres_id_list[[i]]]
    dat$fixed.eff[[i]] <- fixed_eff_mat_sc[ni*(i-1)+pres_id_list[[i]], fixed_eff_ind]
    
    pi_u_i <- exp(dat$fixed.eff[[i]][1,fixed_eff_z_ind] %*% true_params$gamma_u)[1]
    pi_u_i <- pi_u_i/(1+pi_u_i)
    dat$u[i] <- sample(c(1,0), 1, prob=c(pi_u_i, 1-pi_u_i))
    # dat$u[i] <- 0
    pi_w_i <- exp(dat$fixed.eff[[i]][,fixed_eff_w_ind] %*% true_params$gamma_w)
    pi_w_i <- pi_w_i/(1+pi_w_i)
    dat$pi_w[[i]] <- pi_w_i
    flag <- TRUE
    while (flag) {
      dat$w[[i]] <- sapply(
        1:length(pres_id_list[[i]]),
        function(k) { sample(c(1,0), 1, prob=c(pi_w_i[k], 1-pi_w_i[k])) }
      )
      if (sum(dat$w[[i]]) <= 2) { flag <- FALSE }
    }
    dat$w[[i]] <- rep(0, length(pres_id_list[[i]]))
    W_i_mat <- diag(dat$w[[i]])
    I_i <- diag(rep(1, length(pres_id_list[[i]])))
    I_W_i_mat <- I_i - W_i_mat
    pi_z_i <- exp(dat$fixed.eff[[i]][1,fixed_eff_z_ind] %*% true_params$gamma_z)[1]
    pi_z_i <- pi_z_i/(1+pi_z_i)
    dat$z[i] <- sample(c(1,0), 1, prob=c(pi_z_i, 1-pi_z_i))
    # dat$z[i] <- 0
    
    ## sampling sigma_i2
    if (dat$z[i] == 0) {
      dat$sigma_i2[i] <- true_params$sigma_02
      # sigma_i2 <- runif(1, true_params$sigma_02-delta, true_params$sigma_02+delta)
    } else {
      dat$sigma_i2[i] <- runif(
        1,
        true_params$sigma_02*true_params$eta_z^2 - delta*true_params$eta_z^2,
        true_params$sigma_02*true_params$eta_z^2 + delta*true_params$eta_z^2
      )
    }
    
    mean_i <- dat$fixed.eff[[i]][,fixed_eff_ind] %*% true_params$beta
    rand_eff_i <- dat$fixed.eff[[i]][,rand_eff_ind]
    
    flag <- TRUE
    while (flag) {
      r_i <- rmvn_cpp(
        1, 
        rep(0, length(rand_eff_ind)), 
        # true_params$sigma_02 *
        dat$sigma_i2[i] *
          true_params$eta_u^(2*dat$u[i]) * true_params$Lambda,
        1, FALSE
      )[1,]
      eps_i <- rmvn_cpp(
        1, 
        rep(0, length(pres_id_list[[i]])), 
        dat$sigma_i2[i] *
          ((true_params$eta_w-1)*W_i_mat + I_i) %*% 
          h_mat(dat$time[[i]], true_params$rho) %*%
          ((true_params$eta_w-1)*W_i_mat + I_i),
        1, FALSE
      )[1,]
      dat_i <- c(mean_i + rand_eff_i %*% r_i + eps_i)
      # break
      # rand_var_list <- list()
      # eps_var_list <- list()
      if (dat$u[i] == 0 & dat$z[i] == 0 & all(dat$w[[i]] == 0)) {
        flag <- FALSE
      } else {
        # for (k_u in 1:2) {
        #   for (k_z in 1:2) {
        #     rank_var_list[[(k_u-1)*2+k_z]] <- 
        #       true_params$sigma_02 * (true_params$eta_u)^(2*(k_u-1)) *
        #       rand_eff_i %*% true_params$Lambda %*% t(rand_eff_i)
        #     eps_var_list[[(k_u-1)*2+k_z]] <- 
        #       true_params$sigma_02 * true_params$eta_z^(2*(k_z-1)) *
        #       ((true_params$eta_w-1)*W_i_mat + I_i) %*%
        #       h_mat(dat$time[[i]], true_params$rho) %*%
        #       ((true_params$eta_w-1)*W_i_mat + I_i)
        #     llhd <- dmvn_cpp(
        #       matrix(dat_i, nrow = 1),
        #       c(mean_i),
        #       rand_var_list[[(k_u-1)*2+k_z]] + eps_var_list[[(k_u-1)*2+k_z]],
        #       TRUE
        #     )[1,1]
        #   }
        # }
        rand_var_0 <- true_params$sigma_02 *
          rand_eff_i %*% true_params$Lambda %*% t(rand_eff_i)
        eps_var_0 <- true_params$sigma_02 * h_mat(dat$time[[i]], true_params$rho)
        rand_var_i <- true_params$sigma_02 * true_params$eta_z^(2*dat$z[i]) *
          (true_params$eta_u)^(2*dat$u[i]) *
          rand_eff_i %*% true_params$Lambda %*% t(rand_eff_i)
        eps_var_i <- true_params$sigma_02 * true_params$eta_z^(2*dat$z[i]) *
          ((true_params$eta_w-1)*W_i_mat + I_i) %*%
          h_mat(dat$time[[i]], true_params$rho) %*%
          ((true_params$eta_w-1)*W_i_mat + I_i)
        
        diff_llhd <- dmvn_cpp(
          matrix(dat_i, nrow = 1),
          c(mean_i),
          rand_var_i + eps_var_i,
          TRUE
        )[1,1] - dmvn_cpp(
          matrix(dat_i, nrow = 1),
          mean_i,
          rand_var_0 + eps_var_0,
          TRUE
        )[1,1]
        if (diff_llhd > 2) { flag <- FALSE }
      }
    }
    dat$r[[i]] <- r_i
    dat$eps[[i]] <- eps_i
    dat$data[[i]] <- dat_i

    setTxtProgressBar(pb, i)
  }
  dat$r_mat <- do.call(rbind, dat$r)
  return (dat)
}



sim_data_t <- function(
    n, ni, p_miss, true_params, fixed_eff_ind, rand_eff_ind, het_type = "re"
) {
  dat <- list()
  dat$time <- list()
  dat$fixed.eff <- list()
  dat$data <- list()
  dat$r <- list()
  dat$eps <- list()
  pres_id_list <- list()
  for (i in 1:n) {
    miss_id <- sample(c(TRUE,FALSE), ni, replace=TRUE, prob=c(p_miss, 1-p_miss))
    pres_id_list[[i]] <- (1:ni)[!miss_id]
    dat$fixed.eff[[i]] <- cbind(
      1,
      sample(c(0,1),1),
      rnorm(1),
      1:ni,
      (1:ni)^2
    )
  }
  dat$fixed.eff[[i]] <- dat$fixed.eff[[i]][pres_id_list[[i]],]
  fixed_eff_mat <- do.call(rbind, dat$fixed.eff)
  fixed_eff_mat_sc <- cbind(
    1, 
    sapply(
      2:ncol(fixed_eff_mat),
      function(k) { scale(fixed_eff_mat[,k], center=TRUE, scale=TRUE) }
    )
  )
  pb <- txtProgressBar(1,n,style=3)
  for (i in 1:n) {
    dat$time[[i]] <- (1:ni)[pres_id_list[[i]]]
    dat$fixed.eff[[i]] <- fixed_eff_mat_sc[ni*(i-1)+pres_id_list[[i]], fixed_eff_ind] 
    mean_i <- dat$fixed.eff[[i]][,fixed_eff_ind] %*% true_params$beta
    rand_eff_i <- dat$fixed.eff[[i]][,rand_eff_ind]
    
    if (het_type == "re") {
      dat$data[[i]] <- mean_i + rmvt(
        1, 
        true_params$sigma_02*(
          rand_eff_i %*% true_params$Lambda %*% t(rand_eff_i) + 
          h_mat(dat$time[[i]], true_params$rho)
        ), 
        5
      )[1,]
    }
    if (het_type == "r") {
      r_i <- rmvt(1, true_params$sigma_02 * true_params$Lambda, 5)[1,]
      eps_i <- rmvn_cpp(
        1,
        rep(0, pres_id_list[[i]]),
        true_params$sigma_02 * h_mat(dat$time[[i]], true_params$rho),
        1, FALSE
      )[1,]
      dat$data[[i]] <- mean_i + r_i + eps_i
    }
    
    setTxtProgressBar(pb, i)
  }
  return (dat)
}




diagnostic.ar.model <- function(mod, type, burnin, n.sim) {
  ind <- (burnin+1) : (n.sim)
  if (grepl('MI.', type, fixed=TRUE)) {
    if (grepl('component', type, fixed=TRUE)) {
      par(mfrow=c(1,3))
      plot(sqrt(mod$sigma.s2[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[s]))
      plot(sqrt(mod$sigma2[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma))
      plot(mod$rho[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
    } else {
      par(mfrow=c(1,2))
      plot(sqrt(mod$sigma.e2[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma))
      plot(mod$rho[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
    }
  }
  if (grepl('MII.', type, fixed=TRUE)) {
    if (grepl('component', type, fixed=TRUE)) {
      par(mfrow=c(2,3))
      hist(sqrt(mod$sigma.s2[ind,,]), xlab=expression(sigma[s]), breaks=20, main='')
      plot(mod$rho[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
      plot(sqrt(mod$sigma2[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma))
    } else {
      par(mfrow=c(2,3))
      hist(sqrt(mod$sigma.e2[ind,,]), xlab=expression(sigma[e]), breaks=20, main='')
      plot(mod$rho[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
    }
  }
  if (grepl('MIII.', type, fixed=TRUE)) {
    if (grepl('component', type, fixed=TRUE)) {
      par(mfrow=c(2,3))
      hist(sqrt(apply(mod$sigma.s2[ind,,], 2, mean)), xlab=expression(sigma[s]), 
           breaks=20, main='')
      plot(mod$rho[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
      plot(mod$pi.p[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(pi))
      plot(sqrt(mod$sigma2[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma))
    } else {
      par(mfrow=c(2,3))
      hist(sqrt(mod$sigma.e2[ind,,]), xlab=expression(sigma[e]), breaks=20, main='')
      plot(mod$rho[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind,]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
      plot(mod$pi.p[ind,], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(pi))
    }
  }
  par(mfrow=c(1,1))
}

diagnostic.ar.model1 <- function(mod, type, burnin, n.sim, thin) {
  ind <- seq(burnin+1, n.sim, thin)
  if (grepl('MI.', type, fixed=TRUE)) {
    if (grepl('component', type, fixed=TRUE)) {
      par(mfrow=c(1,3))
      plot(sqrt(mod$sigma.s2[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[s]))
      plot(sqrt(mod$sigma2[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma))
      plot(mod$rho[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
    } else {
      par(mfrow=c(1,2))
      plot(sqrt(mod$sigma.e2[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[e]))
      plot(mod$rho[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
    }
  }
  if (grepl('MII.', type, fixed=TRUE)) {
    if (grepl('component', type, fixed=TRUE)) {
      par(mfrow=c(2,3))
      hist(colMeans(sqrt(mod$sigma.s2[ind,])), xlab=expression(sigma[s]), 
           breaks=seq(0,8,0.1), main='')
      plot(mod$rho[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
      plot(sqrt(mod$sigma2[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma))
    } else {
      par(mfrow=c(2,3))
      hist(colMeans(sqrt(mod$sigma.e2[ind,])), xlab=expression(sigma[e]), 
           breaks=seq(0,8,0.1), main='')
      plot(mod$rho[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
    }
  }
  if (grepl('MIII.', type, fixed=TRUE)) {
    if (grepl('component', type, fixed=TRUE)) {
      par(mfrow=c(2,3))
      hist(colMeans(sqrt(mod$sigma.s2[ind,])), xlab=expression(sigma[s]), 
           breaks=seq(0,8,0.1), main='')
      plot(mod$rho[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
      plot(mod$pi.p[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(pi))
      plot(sqrt(mod$sigma2[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma))
    } else {
      par(mfrow=c(2,3))
      hist(colMeans(sqrt(mod$sigma.e2[ind,])), xlab=expression(sigma[e]), 
           breaks=seq(0,20,0.1), main='')
      plot(mod$rho[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
      plot(mod$pi.p[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(pi))
    }
  }
  if (grepl('MIV.', type, fixed=TRUE)) {
    if (grepl('component', type, fixed=TRUE)) {
      par(mfrow=c(3,3))
      hist(colMeans(sqrt(mod$sigma.s2[ind,])), xlab=expression(sigma[s]), 
           breaks=100, main='')
      plot(mod$rho[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(rho))
      plot(mod$alpha[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(alpha))
      plot(sqrt(mod$sigma.02[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[0]))
      plot(sqrt(mod$sigma.12[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma[1]))
      plot(mod$pi.p[ind], type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(pi))
      plot(sqrt(mod$sigma2[ind]), type='l', col='darkgrey', 
           xlab='iteration', ylab=expression(sigma))
    } else {
      plt1 <- ggplot(data.frame(sigma.e2=colMeans(sqrt(mod$sigma.e2[ind,])))) +
        geom_histogram(aes(x=sigma.e2), fill='white', color='darkgrey', bins=100) +
        labs(x=expression(sigma[e])) +
        theme_bw()
      plt2 <- ggplot(data.frame(iter=1:length(mod$rho[ind,]), rho=mod$rho[ind,])) +
        geom_line(aes(x=iter, y=rho), color='darkgrey') +
        labs(y=expression(rho)) +
        theme_bw()
      plt3 <- ggplot(data.frame(iter=1:length(mod$alpha[ind,]), alpha=mod$alpha[ind,])) +
        geom_line(aes(x=iter, y=alpha), color='darkgrey') +
        labs(y=expression(alpha)) +
        theme_bw()
      plt4 <- ggplot(data.frame(iter=1:length(mod$sigma.02[ind,]), 
                                sigma.0=sqrt(mod$sigma.02[ind,]))) +
        geom_line(aes(x=iter, y=sigma.0), color='darkgrey') +
        labs(y=expression(sigma[0])) +
        theme_bw()
      plt5 <- ggplot(data.frame(iter=1:length(mod$sigma.12[ind,]), 
                                sigma.1=sqrt(mod$sigma.12[ind,]))) +
        geom_line(aes(x=iter, y=sigma.1), color='darkgrey') +
        labs(y=expression(sigma[1])) +
        theme_bw()
      plt6 <- ggplot(data.frame(iter=1:length(mod$beta[ind,1]), gam=mod$beta[ind,1])) +
        geom_line(aes(x=iter, y=gam), color='darkgrey') +
        labs(y=expression(beta[0])) +
        theme_bw()
      plt7 <- ggplot(data.frame(iter=1:length(mod$beta[ind,2]), gam=mod$beta[ind,2])) +
        geom_line(aes(x=iter, y=gam), color='darkgrey') +
        labs(y=expression(beta[1])) +
        theme_bw()
      plt8 <- ggplot(data.frame(iter=1:length(mod$beta[ind,3]), gam=mod$beta[ind,3])) +
        geom_line(aes(x=iter, y=gam), color='darkgrey') +
        labs(y=expression(beta[2])) +
        theme_bw()
      plt9 <- ggplot(data.frame(iter=1:length(mod$beta[ind,4]), gam=mod$beta[ind,4])) +
        geom_line(aes(x=iter, y=gam), color='darkgrey') +
        labs(y=expression(beta[3])) +
        theme_bw()
      plt10 <- ggplot(data.frame(iter=1:length(mod$beta[ind,5]), gam=mod$beta[ind,5])) +
        geom_line(aes(x=iter, y=gam), color='darkgrey') +
        labs(y=expression(beta[4])) +
        theme_bw()
      plt11 <- ggplot(data.frame(iter=1:length(mod$gamma[ind,1]), 
                                 gam=mod$gamma[ind,1])) +
        geom_line(aes(x=iter, y=gam), color='darkgrey') +
        labs(y=expression(gamma[0])) +
        theme_bw()
      plt12 <- ggplot(data.frame(iter=1:length(mod$gamma[ind,2]),
                                 gam=mod$gamma[ind,2])) +
        geom_line(aes(x=iter, y=gam), color='darkgrey') +
        labs(y=expression(gamma[3])) +
        theme_bw()
      plt13 <- ggplot(data.frame(iter=1:length(mod$gamma[ind,3]),
                                 gam=mod$gamma[ind,3])) +
        geom_line(aes(x=iter, y=gam), color='darkgrey') +
        labs(y=expression(gamma[4])) +
        theme_bw()
      if (!is.null(mod$pi.e)) {
        plt14 <- ggplot(data.frame(iter=1:length(mod$pi.e[ind,]), pi.e=mod$pi.e[ind,])) +
          geom_line(aes(x=iter, y=pi.e), color='darkgrey') +
          labs(y=expression(pi[e])) +
          theme_bw()
      }
      if (!is.null(mod$pi.e)) {
        print(ggarrange(plt1,plt2,plt3,plt4,plt5,plt6,plt7,plt8,plt9,
                        plt10,plt11,plt12,plt13,plt14, nrow=5, ncol=3))
      } else {
        print(ggarrange(plt1,plt2,plt3,plt4,plt5,plt6,plt7,plt8,plt9, plt10,
                        plt11, plt12, plt13, nrow=5, ncol=3))
      }
      # print(ggarrange(plt1,plt2,plt3,plt4,plt5,plt6,plt7,plt8, plt9,plt10,
      #                 plt11, plt12,plt13, nrow=5, ncol=3))
    }
  }
  if (grepl('MVI.error', type, fixed=TRUE)) {
    plt1 <- ggplot(data.frame(sigma.e2=colMeans(sqrt(mod$sigma.e2[ind,])))) +
      geom_histogram(aes(x=sigma.e2), fill='white', color='darkgrey', bins=100) +
      labs(x=expression(sigma[e])) +
      theme_bw()
    plt2 <- ggplot(data.frame(iter=1:length(mod$rho[ind,]), rho=mod$rho[ind,])) +
      geom_line(aes(x=iter, y=rho), color='darkgrey') +
      labs(y=expression(rho)) +
      theme_bw()
    plt3 <- ggplot(data.frame(iter=1:length(mod$alpha[ind,]), alpha=mod$alpha[ind,])) +
      geom_line(aes(x=iter, y=alpha), color='darkgrey') +
      labs(y=expression(alpha)) +
      theme_bw()
    plt4 <- ggplot(data.frame(iter=1:length(mod$sigma.02[ind,]), 
                              sigma.0=sqrt(mod$sigma.02[ind,]))) +
      geom_line(aes(x=iter, y=sigma.0), color='darkgrey') +
      labs(y=expression(sigma[0])) +
      theme_bw()
    plt5 <- ggplot(data.frame(iter=1:length(mod$sigma.12[ind,]), 
                              sigma.1=sqrt(mod$sigma.12[ind,]))) +
      geom_line(aes(x=iter, y=sigma.1), color='darkgrey') +
      labs(y=expression(sigma[1])) +
      theme_bw()
    plt6 <- ggplot(data.frame(iter=1:length(mod$beta[ind,1]), gam=mod$beta[ind,1])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(beta[0])) +
      theme_bw()
    plt7 <- ggplot(data.frame(iter=1:length(mod$beta[ind,2]), gam=mod$beta[ind,2])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(beta[1])) +
      theme_bw()
    plt8 <- ggplot(data.frame(iter=1:length(mod$beta[ind,3]), gam=mod$beta[ind,3])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(beta[2])) +
      theme_bw()
    plt9 <- ggplot(data.frame(iter=1:length(mod$beta[ind,4]), gam=mod$beta[ind,4])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(beta[3])) +
      theme_bw()
    plt10 <- ggplot(data.frame(iter=1:length(mod$beta[ind,5]), gam=mod$beta[ind,5])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(beta[4])) +
      theme_bw()
    plt11 <- ggplot(data.frame(iter=1:length(mod$gamma.h[ind,1]), 
                               gam=mod$gamma.h[ind,1])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(gamma[h[0]])) +
      theme_bw()
    plt12 <- ggplot(data.frame(iter=1:length(mod$gamma.h[ind,2]),
                               gam=mod$gamma.h[ind,2])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(gamma[h[3]])) +
      theme_bw()
    plt13 <- ggplot(data.frame(iter=1:length(mod$gamma.h[ind,3]),
                               gam=mod$gamma.h[ind,3])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(gamma[h[4]])) +
      theme_bw()
    plt14 <- ggplot(data.frame(iter=1:length(mod$gamma.e[ind,1]), 
                               gam=mod$gamma.e[ind,1])) +
      geom_line(aes(x=iter, y=gam), color='darkgrey') +
      labs(y=expression(gamma[e[0]])) +
      theme_bw()
    plt15 <- 
      ggplot(data.frame(iter=1:length(mod$gamma.e[ind,1]), 
                        pi.e=exp(mod$gamma.e[ind,1])/(1+exp(mod$gamma.e[ind,1])))) +
      geom_line(aes(x=iter, y=pi.e), color='darkgrey') +
      labs(y=expression(pi[e])) +
      theme_bw()
    # not printing this plot
    plt16 <- ggplot(data.frame(pi.h=colMeans(mod$pi.h[ind,]))) +
      geom_histogram(aes(x=pi.h), fill='white', color='darkgrey', bins=100) +
      labs(x=expression(pi[h])) +
      theme_bw()
    print(ggarrange(plt1,plt2,plt3,plt4,plt5,plt6,plt7,plt8, 
                    plt9,plt10,
                    plt11, 
                    plt12,plt13,
                    plt14,plt15,
                    nrow=5, ncol=3))
  }
  par(mfrow=c(1,1))
}

plot.sigma <- function(mcmc.df, ind) {
  cls <- mcmc.df$cls[ind,]
  sigma.e <- sqrt(mcmc.df$sigma.e2[ind,])
  K <- ncol(sigma.e)
  for (i in 1:length(ind)) {
    aux.cls.ids <- which(!(1:K %in% unlist(cls[i,])))
    sigma.e[i, aux.cls.ids] <- NA
  }
  plot.df <- data.frame(iter = rep(1:length(ind), times=K),
                        sigma.e = c(data.matrix(sigma.e)),
                        cluster = factor(rep(1:K, each=length(ind))))
  plt <- ggplot(plot.df, aes(x=iter, y=sigma.e)) +
    geom_point(aes(color=cluster), size=0) +
    geom_line(aes(color=cluster), alpha=0.6) +
    labs(y=expression(sigma[e])) +
    guides(color="none") +
    theme_bw()
  return(list(ggMarginal(plt, plot.df, type="density", margins="y", 
                         fill = alpha("salmon", 0.3))))
}

plt.mcmc <- function(df, n.chain, param.ind, param.type, param.label, have.subscript,
                     param.label.x = NULL, x.range = NULL, y.range = NULL, 
                     title = NULL) {
  if (n.chain == 1) {
    if (is.vector(df)) df <- as.matrix(df)
    if (ncol(df) == 1 & !is.null(param.ind))
      stop("error: param.ind must be NULL if df is a vector or column matrix.\n")
  }
  param.len <- ifelse (n.chain == 1, ncol(df), ncol(df)/n.chain)
  if (is.null(param.ind)) {
    param.ind <- 1:param.len
  }
  if (param.type == "histogram") {
    combined.param.ind <- c()
    for (i.chain in 1:n.chain) {
      combined.param.ind <- c(combined.param.ind, (i.chain-1)*param.len + param.ind)
    }
    hist.df <- data.frame(x = colMeans(df)[combined.param.ind],
                          type = factor(rep(1:n.chain, each = length(param.ind))))
    plt <- 
      ggplot(hist.df, aes(x = x, fill = type)) +
      geom_density(aes(y = ..density..), alpha = 0.25) +
      labs(x=parse(text = param.label)) +
      guides(fill = "none") +
      theme_bw()
    plt.list <- list(plt)
  } else {
    plt.list <- list()
    for (j in param.ind) {
      combined.param.ind.j <- c()
      for (i.chain in 1:n.chain) {
        combined.param.ind.j <- c(combined.param.ind.j, (i.chain-1)*param.len + j)
      }
      chain.df <- data.frame(y = c(df[,combined.param.ind.j]),
                             x = rep(1:nrow(df), times = n.chain),
                             type = factor(rep(1:n.chain, each = nrow(df))))
      plt <- ggplot(chain.df, aes(x=x, y=y, color=type))+
        geom_point(size = 0) +
        geom_line(alpha = 0.6) +
        guides(color = "none") +
        theme_bw()
      if (length(param.ind) == 1 & !have.subscript) {
        plt <- plt + labs(x='iter', y=parse(text = param.label))
      } else {
        y_lab <- ifelse (is.null(param.label.x),
                         paste0(param.label, '[', j-1, ']'),
                         paste0(param.label, '[', j-1, ']', 
                                ' (', param.label.x[j], ')'))
        plt <- plt + labs(x='iter', 
                          y=parse(text = y_lab))
      }
      plt.list[[j]] <- plt
    }
  }
  
  if (!is.null(x.range)) {
    for (i in 1:length(plt.list)) {
      plt.list[[i]] <- plt.list[[i]] + lims(x = x.range) 
    }
  }
  if (!is.null(y.range)) {
    for (i in 1:length(plt.list)) {
      plt.list[[i]] <- plt.list[[i]] + lims(y = y.range) 
    }
  }
  if (!is.null(title)) {
    for (i in 1:length(plt.list)) {
      plt.list[[i]] <- 
        plt.list[[i]] + 
        labs(tag = parse(text = title)) +
        theme(plot.tag.position = 'bottom')
      # theme(plot.title = element_text(hjust=0.5, size=10))
    }
  }
  return (plt.list)
}

plt.model <- function(mod.list, n.chain, ind, type, sigma.range=NULL,
                      beta.labels = NULL, gam.e.labels = NULL, gam.h.labels = NULL,
                      gam.W.labels = NULL) {
  param.name.vec <- names(mod.list[[1]])
  n.param <- length(param.name.vec)
  mod <- list()
  for (param.name in param.name.vec) {
    mod[[param.name]] <- array()
    for (i.chain in 1:n.chain) {
      mod[[param.name]] <- cbind(mod[[param.name]], mod.list[[i.chain]][[param.name]])
    }
    mod[[param.name]] <- mod[[param.name]][,-1, drop=FALSE] ## removing NA column
  }
  plt.list <- list()
  if (type == '00') {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.12[ind,]), n.chain, NULL, '', 'sigma[1]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == '01') {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == '10') {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.12[ind,]), n.chain, NULL, '', 'sigma[1]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == '11') {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == 'IV') {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(sqrt(mod$sigma.12[ind,]), n.chain, NULL, '', 'sigma[1]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(mod$pi.h[ind,], n.chain, NULL, 'histogram', 'pi[h]', 
                           FALSE),
                  plt.mcmc(mod$gamma[ind,], n.chain, NULL, '', 'gamma[h]', 
                           TRUE, gam.h.labels),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == 'IV_II') {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(mod$pi.h[ind,], n.chain, NULL, 'histogram', 'pi[h]', 
                           FALSE),
                  plt.mcmc(mod$gamma[ind,], n.chain, NULL, '', 'gamma[h]', 
                           TRUE, gam.h.labels),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == 'V') {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(sqrt(mod$sigma.12[ind,]), n.chain, NULL, '', 'sigma[1]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(mod$pi.e[ind,], n.chain, NULL, '', 'pi[e]', FALSE),
                  plt.mcmc(mod$pi.h[ind,], n.chain, NULL, 'histogram', 'pi[h]', FALSE),
                  plt.mcmc(mod$gamma[ind,], n.chain, NULL, '', 'gamma[h]', 
                           TRUE, gam.h.labels),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == 'VI') {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(sqrt(mod$sigma.12[ind,]), n.chain, NULL, '', 'sigma[1]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(mod$pi.e[ind,], n.chain, NULL, 'histogram', 'pi[e]', FALSE),
                  plt.mcmc(mod$pi.h[ind,], n.chain, NULL, 'histogram', 'pi[h]', FALSE),
                  plt.mcmc(mod$gamma.e[ind,], n.chain, NULL, '', 'gamma[e]', 
                           TRUE, gam.e.labels),
                  plt.mcmc(mod$gamma.h[ind,], n.chain, NULL, '', 'gamma[h]', 
                           TRUE, gam.h.labels),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == "DPM") {
    plt1 <- plot.sigma(mod, ind)
    plt.list <- c(plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
    n.plt <- length(plt.list)
    plt.grid <- matrix(NA, nrow=ceiling(n.plt/3), ncol=3)
    plt.grid[1:n.plt] <- 2:(n.plt+1)
    print(grid.arrange(grobs = as.list(c(plt1, plt.list)),
                       width = c(1,1,1),
                       layout_matrix = rbind(c(1,1,1),
                                             c(1,1,1),
                                             t(plt.grid))))
    return ()
  } else if (type == "VII") {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(sqrt(mod$sigma.12[ind,]), n.chain, NULL, '', 'sigma[1]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(sqrt(mod$sigma.22[ind,]), n.chain, NULL, '', 'sigma[2]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(mod$pi.h[ind,], n.chain, NULL, 'histogram', 'pi[h]', 
                           FALSE),
                  plt.mcmc(mod$pi.W[ind,], n.chain, NULL, 'histogram', 'pi[W]', 
                           FALSE),
                  plt.mcmc(mod$gamma.h[ind,], n.chain, NULL, '', 'gamma[h]', 
                           TRUE, gam.h.labels),
                  plt.mcmc(mod$gamma.W[ind,], n.chain, NULL, '', 'gamma[W]', 
                           TRUE, gam.W.labels),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == "VII_I") {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(sqrt(mod$sigma.22[ind,]), n.chain, NULL, '', 'sigma[2]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  # plt.mcmc(mod$pi.h[ind,], n.chain, NULL, 'histogram', 'pi[h]', 
                  #          FALSE),
                  plt.mcmc(rowMeans(mod$z[ind,]), n.chain, NULL, '', 'bar(z)', FALSE),
                  # plt.mcmc(mod$pi.W[ind,], n.chain, NULL, 'histogram', 'pi[W]', 
                  #          FALSE),
                  plt.mcmc(rowMeans(mod$W[ind,], na.rm=TRUE), n.chain, NULL, '', 
                           'bar(W)', FALSE),
                  plt.mcmc(mod$gamma.h[ind,], n.chain, NULL, '', 'gamma[h]', 
                           TRUE, gam.h.labels),
                  plt.mcmc(mod$gamma.W[ind,], n.chain, NULL, '', 'gamma[W]', 
                           TRUE, gam.W.labels),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  } else if (type == "VIII") {
    plt.list <- c(plt.mcmc(sqrt(mod$sigma.e2[ind,]), n.chain, NULL, 'histogram',
                           'sigma[e]', FALSE),
                  plt.mcmc(mod$alpha[ind,], n.chain, NULL, '', 'alpha', FALSE),
                  plt.mcmc(sqrt(mod$sigma.02[ind,]), n.chain, NULL, '', 'sigma[0]', FALSE,
                           y.range = sigma.range),
                  plt.mcmc(mod$rho[ind,], n.chain, NULL, '', 'rho', FALSE),
                  plt.mcmc(rowMeans(mod$z[ind,]), n.chain, NULL, '', 'bar(z)', FALSE),
                  plt.mcmc(rowMeans(mod$W[ind,], na.rm=TRUE), n.chain, NULL, '', 
                           'bar(W)', FALSE),
                  plt.mcmc(mod$gamma.h[ind,], n.chain, NULL, '', 'gamma[h]', 
                           TRUE, gam.h.labels),
                  plt.mcmc(mod$gamma.W[ind,], n.chain, NULL, '', 'gamma[W]', 
                           TRUE, gam.W.labels),
                  plt.mcmc(mod$beta[ind,], n.chain, NULL, '', 'beta', TRUE, beta.labels))
  }
  
  print(ggarrange(plotlist = plt.list, 
                  nrow=ceiling(length(plt.list)/3), 
                  ncol=3))
}

diagnostic_dpm <- function(dat, mod, ind) {
  # ggplot(data.frame(x = 1:length(ind),
  #                   y = sapply(ind, 
  #                              function(iter) {
  #                                sqrt(mod$sigma.e2[iter, mod$cls[iter, i]])
  #                              }))) +
  #   geom_line(aes(x=x, y=y), color="steelblue") +
  #   theme_bw()
  mean_pred_mat <- matrix(NA, nrow = nrow(dat$data), ncol = ncol(dat$data))
  for (i in 1:nrow(dat$data)) {
    mean_pred_mat[i,] <- colMeans(t(sapply(ind,
                                           function(iter) {
                                             dat$fixed.eff[[i]] %*% mod$beta[iter,]
                                           })))
  }
  return(ggplot(data.frame(data = c(dat$data),
                           fitted_data = c(mean_pred_mat))) +
           geom_point(aes(x = data, y = fitted_data), color="steelblue", size=0.5) +
           geom_abline(slope = 1, intercept = 0, color="salmon", size=1) +
           labs(x = expression(Y[ij]), y = expression(hat(Y)[ij])) +
           theme_bw())
}

predictive.check <- function(dat) {
  n <- nrow(dat)
  ni.max <- ncol(dat)
  profile.mean <- rowMeans(dat)
  temp.dat <- dat - matrix(rep(profile.mean, times=ni.max), nrow=n, byrow = FALSE)
  paired.dat <- matrix(NA, nrow=n, ncol=ni.max-1)
  for (j in 1:(ni.max-1)) {
    paired.dat[,j] <- temp.dat[,j] * temp.dat[,j+1]
  }
  pairwise.cov.vec <- rowSums(paired.dat)/(ni.max-1)
  var1.dat <- var2.dat <- matrix(NA, nrow=n, ncol=ni.max-1)
  for (j in 1:(ni.max-1)) {
    var1.dat[,j] <- temp.dat[,j]^2
    var2.dat[,j] <- temp.dat[,(j+1)]^2
  }
  var.dat <- cbind(var1.dat, var2.dat[,(ni.max-1)])
  
  var1.vec <- rowSums(var1.dat) / (ni.max - 1)
  var2.vec <- rowSums(var2.dat) / (ni.max - 1)
  var.vec <- rowSums(var.dat) / (ni.max - 1)
  # corr.vec <- pairwise.cov.vec / sqrt(var1.vec * var2.vec)
  corr.vec <- pairwise.cov.vec / var.vec
  return (list(var.stat = var.vec,
               cov.stat = pairwise.cov.vec,
               corr.stat = corr.vec))
}

predictive.check <- function(dat) {
  n <- length(dat)
  temp.dat <- list()
  pairwise.cov.vec <- var.vec <- array(NA, n)
  for (i in 1:n) {
    temp.dat[[i]] <- dat[[i]]
    pairwise.cov.vec[i] <- 
      sum(temp.dat[[i]][-1] * temp.dat[[i]][-length(temp.dat[[i]])]) / 
      (length(temp.dat[[i]]) - 1)
    var.vec[i] <-
      sum(temp.dat[[i]]^2) / (length(temp.dat[[i]] - 1))
  }
  corr.vec <- pairwise.cov.vec / var.vec
  return (list(var.stat = var.vec,
               cov.stat = pairwise.cov.vec,
               corr.stat = corr.vec))
}

predictive_diag <- function(dat, model, model.type) {
  ind <- seq(5001, 10000, 5)
  
  # n <- nrow(dat$data)
  n <- length(dat$data)
  # ni_max <- ncol(dat$data)
  
  corr.stat <- corr.stat_centered <- array(NA, length(ind))
  true.corr.stat <- array(NA, length(ind))
  
  # cov.stat <- cov.stat_centered <- array(NA, dim=c(length(ind), nrow(dat$data)))
  # var.stat <- var.stat_centered <- array(NA, dim=c(length(ind), nrow(dat$data)))
  cov.stat <- cov.stat_centered <- array(NA, dim=c(length(ind), length(dat$data)))
  var.stat <- var.stat_centered <- array(NA, dim=c(length(ind), length(dat$data)))
  
  true.corr.stat_centered <- array(NA, length(ind))
  
  # true.cov.stat_centered <- array(NA, dim=c(length(ind), nrow(dat$data)))
  # true.var.stat_centered <- array(NA, dim=c(length(ind), nrow(dat$data)))
  true.cov.stat_centered <- array(NA, dim=c(length(ind), length(dat$data)))
  true.var.stat_centered <- array(NA, dim=c(length(ind), length(dat$data)))
  
  # true_data_centered <- matrix(NA, nrow = n, ncol = ni_max)
  true_data_centered <- list()
  
  for (iter in 1:length(ind)) {
    if (model.type == "DPM") {
      beta.s <- model$beta[ind[iter],]
      K <- ncol(model$sigma.e2)
      cls.s <- model$cls[ind[iter],]
      prob.cls.s <- model$prob.cls[ind[iter],]
      sigma.e2.s <- model$sigma.e2[ind[iter],]
      rho.s <- model$rho[ind[iter]]
    } else if (model.type == "IV" | model.type == "IV_I") {
      beta.s <- model$beta[ind[iter],]
      sigma.02.s <- model$sigma.02[ind[iter]]
      sigma.12.s <- model$sigma.12[ind[iter]]
      pi.h.s <- model$pi.h[ind[iter],]
      alpha.s <- model$alpha[ind[iter]]
      sigma.e2.s <- model$sigma.e2[ind[iter],]
      rho.s <- model$rho[ind[iter]]
    } else if (model.type == "00") {
      beta.s <- model$beta[ind[iter],]
      sigma.12.s <- model$sigma.12[ind[iter]]
      alpha.s <- model$alpha[ind[iter]]
      sigma.e2.s <- model$sigma.e2[ind[iter],]
      rho.s <- model$rho[ind[iter]]
    }
    
    # data_sim <- data_sim_centered <- matrix(NA, nrow = n, ncol = ni_max)
    data_sim <- data_sim_centered <- list()
    
    # for (i in 1:nrow(dat$data)) {
    for (i in 1:length(dat$data)) {
      tt <- dat$time[[i]]
      if (model.type == "DPM") {
        cls.s[i] <- sample(1:K, 1, prob = prob.cls.s[((i-1)*K+1):(i*K)])
        # data_sim[i, tt] <- mvrnorm(1,
        #                            dat$fixed.eff[[i]][tt,] %*% beta.s,
        #                            sigma.e2.s[cls.s[i]] * h_mat(dat$time[[i]], rho.s))
        data_sim[[i]] <- mvrnorm(1,
                                 dat$fixed.eff[[i]] %*% beta.s,
                                 sigma.e2.s[cls.s[i]] * h_mat(dat$time[[i]], rho.s))
      } else if (model.type == "IV") {
        z <- sample(c(1,0), 1, prob = c(pi.h.s[i], 1-pi.h.s[i]))
        sigma.e2.s[i] <- ifelse (z == 1,
                                 sigma.02.s,
                                 rgamma(1, alpha.s, alpha.s/sigma.12.s))
        # data_sim[i, tt] <- mvrnorm(1,
        #                            dat$fixed.eff[[i]][tt,] %*% beta.s,
        #                            sigma.e2.s[i] * h_mat(dat$time[[i]], rho.s))
        data_sim[[i]] <- mvrnorm(1,
                                 dat$fixed.eff[[i]] %*% beta.s,
                                 sigma.e2.s[i] * h_mat(dat$time[[i]], rho.s))
      } else if (model.type == "IV_I") {
        z <- sample(c(1,0), 1, prob = c(pi.h.s[i], 1-pi.h.s[i]))
        sigma.e2.s[i] <- ifelse (z == 1,
                                 sigma.02.s,
                                 rgamma(1, alpha.s+1, alpha.s/sigma.12.s))
        # data_sim[i, tt] <- mvrnorm(1,
        #                            dat$fixed.eff[[i]][tt,] %*% beta.s,
        #                            sigma.e2.s[i] * h_mat(dat$time[[i]], rho.s))
        data_sim[[i]] <- mvrnorm(1,
                                 dat$fixed.eff[[i]] %*% beta.s,
                                 sigma.e2.s[i] * h_mat(dat$time[[i]], rho.s))
      } else if (model.type == "00") {
        sigma.e2.s[i] <- rgamma(1, alpha.s, alpha.s/sigma.12.s)
        data_sim[[i]] <- mvrnorm(1,
                                 dat$fixed.eff[[i]] %*% beta.s,
                                 sigma.e2.s[i] * h_mat(dat$time[[i]], rho.s))
      }
    }
    # for (i in 1:nrow(dat$data)) {
    for (i in 1:length(dat$data)) {
      tt <- dat$time[[i]]
      # data_sim_centered[i, tt] <- 
      #   data_sim[i, tt] - dat$fixed.eff[[i]][tt,] %*% beta.s
      # true_data_centered[i, tt] <-
      #   dat$data[i, tt] - dat$fixed.eff[[i]][tt,] %*% beta.s
      data_sim_centered[[i]] <- 
        data_sim[[i]] - dat$fixed.eff[[i]] %*% beta.s
      true_data_centered[[i]] <-
        dat$data[[i]] - dat$fixed.eff[[i]] %*% beta.s
    }
    
    res <- predictive.check(data_sim)
    res_centered <- predictive.check(data_sim_centered)
    
    corr.stat[iter] <- mean(res$corr.stat)
    corr.stat_centered[iter] <- mean(res_centered$corr.stat)
    
    var.stat[iter,] <- res$var.stat
    var.stat_centered[iter,] <- res_centered$var.stat
    
    cov.stat[iter,] <- res$cov.stat
    cov.stat_centered[iter,] <- res_centered$cov.stat
    
    true_res_centered <- predictive.check(true_data_centered)
    true.corr.stat_centered[iter] <- mean(true_res_centered$corr.stat)
    true.cov.stat_centered[iter,] <- true_res_centered$cov.stat
    true.var.stat_centered[iter,] <- true_res_centered$var.stat
  }
  
  return (list(corr.stat = corr.stat,
               var.stat = var.stat,
               cov.stat = cov.stat,
               true.corr.stat = predictive.check(dat$data)$corr.stat,
               true.var.stat = predictive.check(dat$data)$var.stat,
               true.cov.stat = predictive.check(dat$data)$cov.stat,
               corr.stat_centered = corr.stat_centered,
               var.stat_centered = var.stat_centered,
               cov.stat_centered = cov.stat_centered,
               true.corr.stat_centered = true.corr.stat_centered,
               true.var.stat_centered = true.var.stat_centered,
               true.cov.stat_centered = true.cov.stat_centered))
}

plt.predictive.checks <- function(diag.list) {
  # hist(diag.list$corr.stat - mean(diag.list$true.corr.stat), 
  #      breaks=100, main="Correlation: Estimated - True",
  #      xlab = "Difference in Correlation")
  plt.cor <- 
    ggplot(data.frame(x = diag.list$corr.stat_centered - 
                        diag.list$true.corr.stat_centered)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Correlation",
         title = "Correlation: Estimated - True (Centered)") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.05 <- 
    apply(sqrt(diag.list$var.stat_centered), 1, quantile, 0.05)
  true.sd.stat_centered.05 <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, quantile, 0.05)
  plt.sd.05 <- 
    ggplot(data.frame(x = sd.stat_centered.05 - true.sd.stat_centered.05)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n 5% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.20 <- 
    apply(sqrt(diag.list$var.stat_centered), 1, quantile, 0.20)
  true.sd.stat_centered.20 <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, quantile, 0.20)
  plt.sd.20 <- 
    ggplot(data.frame(x = sd.stat_centered.20 - true.sd.stat_centered.20)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n 20% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.30 <- 
    apply(sqrt(diag.list$var.stat_centered), 1, quantile, 0.30)
  true.sd.stat_centered.30 <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, quantile, 0.30)
  plt.sd.30 <- 
    ggplot(data.frame(x = sd.stat_centered.30 - true.sd.stat_centered.30)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n 30% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.40 <- 
    apply(sqrt(diag.list$var.stat_centered), 1, quantile, 0.40)
  true.sd.stat_centered.40 <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, quantile, 0.40)
  plt.sd.40 <- 
    ggplot(data.frame(x = sd.stat_centered.40 - true.sd.stat_centered.40)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n 40% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.median <- 
    apply(sqrt(diag.list$var.stat_centered), 1, median)
  true.sd.stat_centered.median <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, median)
  plt.sd.med <- 
    ggplot(data.frame(x = sd.stat_centered.median - true.sd.stat_centered.median)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n Median") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.60 <- 
    apply(sqrt(diag.list$var.stat_centered), 1, quantile, 0.60)
  true.sd.stat_centered.60 <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, quantile, 0.60)
  plt.sd.60 <- 
    ggplot(data.frame(x = sd.stat_centered.60 - true.sd.stat_centered.60)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n 60% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.70 <- 
    apply(sqrt(diag.list$var.stat_centered), 1, quantile, 0.70)
  true.sd.stat_centered.70 <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, quantile, 0.70)
  plt.sd.70 <- 
    ggplot(data.frame(x = sd.stat_centered.70 - true.sd.stat_centered.70)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n 70% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.80 <- 
    apply(sqrt(diag.list$var.stat_centered), 1, quantile, 0.80)
  true.sd.stat_centered.80 <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, quantile, 0.80)
  plt.sd.80 <- 
    ggplot(data.frame(x = sd.stat_centered.80 - true.sd.stat_centered.80)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n 80% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  sd.stat_centered.95 <- 
    apply(sqrt(diag.list$var.stat_centered), 1, quantile, 0.95)
  true.sd.stat_centered.95 <- 
    apply(sqrt(diag.list$true.var.stat_centered), 1, quantile, 0.95)
  plt.sd.95 <- 
    ggplot(data.frame(x = sd.stat_centered.95 - true.sd.stat_centered.95)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Std. Dev.",
         title = "Std. Dev.: Estimated - True (Centered) \n 95% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.05 <- 
    apply(diag.list$cov.stat_centered, 1, quantile, 0.05)
  true.cov.stat_centered.05 <- 
    apply(diag.list$true.cov.stat_centered, 1, quantile, 0.05)
  plt.cov.05 <- 
    ggplot(data.frame(x = cov.stat_centered.05 - true.cov.stat_centered.05)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n 5% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.20 <- 
    apply(diag.list$cov.stat_centered, 1, quantile, 0.20)
  true.cov.stat_centered.20 <- 
    apply(diag.list$true.cov.stat_centered, 1, quantile, 0.20)
  plt.cov.20 <- 
    ggplot(data.frame(x = cov.stat_centered.20 - true.cov.stat_centered.20)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n 20% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.30 <- 
    apply(diag.list$cov.stat_centered, 1, quantile, 0.30)
  true.cov.stat_centered.30 <- 
    apply(diag.list$true.cov.stat_centered, 1, quantile, 0.30)
  plt.cov.30 <- 
    ggplot(data.frame(x = cov.stat_centered.30 - true.cov.stat_centered.30)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n 30% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.40 <- 
    apply(diag.list$cov.stat_centered, 1, quantile, 0.40)
  true.cov.stat_centered.40 <- 
    apply(diag.list$true.cov.stat_centered, 1, quantile, 0.40)
  plt.cov.40 <- 
    ggplot(data.frame(x = cov.stat_centered.40 - true.cov.stat_centered.40)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n 40% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.median <- 
    apply(diag.list$cov.stat_centered, 1, median)
  true.cov.stat_centered.median <- 
    apply(diag.list$true.cov.stat_centered, 1, median)
  plt.cov.med <- 
    ggplot(data.frame(x = cov.stat_centered.median - true.cov.stat_centered.median)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n Median") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.60 <- 
    apply(diag.list$cov.stat_centered, 1, quantile, 0.60)
  true.cov.stat_centered.60 <- 
    apply(diag.list$true.cov.stat_centered, 1, quantile, 0.60)
  plt.cov.60 <- 
    ggplot(data.frame(x = cov.stat_centered.60 - true.cov.stat_centered.60)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n 60% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.70 <- 
    apply(diag.list$cov.stat_centered, 1, quantile, 0.70)
  true.cov.stat_centered.70 <- 
    apply(diag.list$true.cov.stat_centered, 1, quantile, 0.70)
  plt.cov.70 <- 
    ggplot(data.frame(x = cov.stat_centered.70 - true.cov.stat_centered.70)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n 70% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.80 <- 
    apply(diag.list$cov.stat_centered, 1, quantile, 0.80)
  true.cov.stat_centered.80 <- 
    apply(diag.list$true.cov.stat_centered, 1, quantile, 0.80)
  plt.cov.80 <- 
    ggplot(data.frame(x = cov.stat_centered.80 - true.cov.stat_centered.80)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n 80% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  cov.stat_centered.95 <- 
    apply(diag.list$cov.stat_centered, 1, quantile, 0.95)
  true.cov.stat_centered.95 <- 
    apply(diag.list$true.cov.stat_centered, 1, quantile, 0.95)
  plt.cov.95 <- 
    ggplot(data.frame(x = cov.stat_centered.95 - true.cov.stat_centered.95)) +
    geom_histogram(aes(x=x), bins=50, color='black', fill='white') +
    geom_vline(xintercept = 0, color="red") +
    labs(x="(Estimated - True) Covariance",
         title = "Covariance: Estimated - True (Centered) \n 95% Quantile") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))
  
  ggarrange(plotlist = list(plt.cor, plt.sd.05, plt.sd.20, plt.sd.30, plt.sd.40, 
                            plt.sd.med, plt.sd.60, plt.sd.70, plt.sd.80, plt.sd.95,
                            plt.cov.05, plt.cov.20, plt.cov.30, plt.cov.40, plt.cov.med, 
                            plt.cov.60,  plt.cov.70, plt.cov.80, plt.cov.95),
            nrow = 7, ncol=3)
}

diagnostic.data <- function(dat) {
  n = nrow(dat$data); ni = ncol(dat$data)
  if (is.null(dat$time)) { tt <- dat$time } else { tt <- 1:ni }
  data.df <- data.frame(Subject = rep(1:n, each=ni), 
                        Time = rep(tt, times=n), 
                        Y = c(t(dat$data)))
  profile.plot <- 
    ggplot(data.df) + 
    geom_point(aes(x=Time, y=Y, group=Subject), color="black", size=0.8) +
    geom_line(aes(x=Time, y=Y, group=Subject), color="steelblue") +
    scale_x_continuous(name="Age (in years)", breaks=1:ni, minor_breaks = NULL) +
    theme_bw()
  profile.density.plot <-
    ggMarginal(profile.plot, data.df, type="density", margins="y", color="black", 
               fill="steelblue")
  rand.eff.df <- data.frame(cbreaksd(1:n, dat$rand.eff))
  names(rand.eff.df) <- c("Subject", "b.0", "b.1")
  hist.rand.eff.plot <-
    ggplot(rand.eff.df) + geom_point(aes(x=b.0, y=b.1), color="black", size=1) +
    scale_x_continuous(name="Random intercept") +
    scale_y_continuous(name="Time effect") +
    theme_bw()
  rand.eff.plot <- 
    ggMarginal(hist.rand.eff.plot, rand.eff.df, type="density", margins="both", 
               color="black", fill="steelblue")
  serial.eff.df <- data.frame(c(t(dat$serial.eff))); names(serial.eff.df) <- c("w")
  hist.serial.eff.plot <-
    ggplot(serial.eff.df) + geom_histogram(aes(x=w), fill="steelblue", color="black", 
                                           breakss=50) +
    scale_x_continuous(name="Serial correlation") +
    theme_bw()
  profile.density.plot <- annotate_figure(profile.density.plot, bottom="(a)", 
                                          fig.lab.pos="bottom")
  rand.eff.plot <- annotate_figure(rand.eff.plot, bottom="(b)", fig.lab.pos="top")
  hist.serial.eff.plot <- annotate_figure(hist.serial.eff.plot, bottom = "(c)", 
                                          fig.lab.pos = "top")
  print(ggarrange(profile.density.plot, rand.eff.plot, hist.serial.eff.plot, ncol=3))
}

compute.LPML <- function(model, dat, n.sim, n.burn, thin) {
  ind <- seq(n.burn+1, n.sim, thin)
  LPML <- 0
  for (i in 1:nrow(dat$data)) {
    cpo.i <- 0
    for (j in ind) {
      mean.y <- dat$fixed.eff[[i]][dat$time[[i]],] %*% model$beta[j,]
      Sigma.y <- H.mat(length(dat$time[[i]]), 
                       dat$time[[i]], 
                       sqrt(model$sigma.e2[j,i]),
                       model$rho[j])
      cpo.i <- cpo.i + 1/dmvnorm(dat$data[i,dat$time[[i]]], mean.y, Sigma.y, log=FALSE)
    }
    cpo.i <- cpo.i/length(ind)
    LPML <- LPML - log(cpo.i)
  }
  return (LPML)
}



compute_ll <- function(
    dat, mod, iter_arr,
    fixed_eff_ind, rand_eff_ind,
    fixed_eff_u_ind, fixed_eff_w_ind, fixed_eff_z_ind
) {
  ni <- sapply(dat$time, length)
  n <- length(ni)
  ni_max <- max(ni)
  ll_mat <- matrix(NA, nrow = length(iter_arr), ncol = n)
  pb <- txtProgressBar(1, length(iter_arr), style=3)
  for (i_iter in seq_along(iter_arr)) {
    iter <- iter_arr[i_iter]
    beta <- mod$mcmc.df$beta[iter,]
    gamma.u <- mod$mcmc.df$gamma.u[iter,]
    gamma.z <- mod$mcmc.df$gamma.h[iter,]
    gamma.w <- mod$mcmc.df$gamma.W[iter,]
    alpha_arr <- mod$mcmc.df$alpha[iter,]
    sigma_02 <- mod$mcmc.df$sigma.02[iter]
    rho <- mod$mcmc.df$rho[iter]
    Lambda <- matrix(mod$mcmc.df$Lambda[iter,], nrow=2)
    for (i in 1:n) {
      fixed_eff_i <- dat$fixed.eff[[i]][,fixed_eff_ind]
      rand_eff_i <- dat$fixed.eff[[i]][,rand_eff_ind]
      fixed_eff_u <- dat$fixed.eff[[i]][,fixed_eff_u_ind]
      fixed_eff_w <- dat$fixed.eff[[i]][,fixed_eff_w_ind]
      fixed_eff_z <- dat$fixed.eff[[i]][,fixed_eff_z_ind]
      u_i <- mod$mcmc.df$u[i]
      z_i <- mod$mcmc.df$z[i]
      W_i <- mod$mcmc.df$W[iter, (ni_max*(i-1)+1) : (ni_max*i)]
      W_i_mat <- diag(W_i)
      I_W_i_mat <- diag(1, ni_max) - W_i_mat
      H_inv <- solve(
        (eta_u)^(2*u_i) * (rand_eff_i %*% Lambda %*% t(rand_eff_i)) +
          (eta_w*W_i_mat + I_W_i_mat) %*% 
          h_mat(dat$time[[i]], rho) %*% 
          (eta_w*W_i_mat + I_W_i_mat)
      )
      quad <- 
        c(t(dat$data[[i]] - fixed_eff_i %*% beta) %*%
            H_inv %*%
            (dat$data[[i]] - fixed_eff_i %*% beta))
      zeta <- 1/alpha_arr[z_i+1]^2 - ni[i]/2
      ll_mat[i_iter, i] <- compute_pdf(
        quad, alpha_arr, z_i, eta_z, sigma_02, H_inv, ni[i]
      )
    }
    setTxtProgressBar(pb, i_iter)
  }
  return (ll_mat)
}


