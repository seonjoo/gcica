#' Implement different W same T group-colorICA algorithm.
#'
#' param Xc      Hierachical list containing       [list]
#' data matrices. G is the number of groups and NG is the number of subjects in each group.
#' Each component of the list is a list with data matrices with dimension M by N,
#' where M is the number of mixtures and N is the number of time points.
#' Each data matrix has been pre-whitenned and centered.
#'
#' return W: a list of M by M unmixing matries
#' export
#'
#' examples
#' 
#################
# Load packages #
#################
.libPaths('/ifs/scratch/msph/LeeLab/software/R/hpc') # set the R library
library(parallel)
library(pracma)
library(itsmr)
library(coloredICA)
n_core = 10

#################################################
# Modified Yule-Walker Algorithm for Group Data #
#################################################


group_yw = function(Xc, maxnmodels = 10){
  all_models = mclapply(1:maxnmodels, function(p){
    n = length(Xc[1,])
    gamma = lapply(1:nrow(Xc),function(i){
      x_center = Xc[i,] - mean(Xc[i,])
      gamma = acvf(x_center, p)
      return(gamma)
    })
    Gamma = lapply(1:nrow(Xc),function(i){
      x_center = Xc[i,] - mean(Xc[i,])
      gamma = acvf(x_center, p)
      Gamma = toeplitz(gamma[1:p])
      return(Gamma)
    })
    gamma = Reduce("+", gamma) / length(gamma)
    Gamma = Reduce("+", Gamma) / length(Gamma)
    phi = solve(Gamma, gamma[2:(p + 1)])
    v = gamma[1] - drop(crossprod(gamma[2:(p + 1)], phi))
    V = v * solve(Gamma)
    se.phi = sqrt(1/n * diag(V))
    a = list(phi = phi, theta = 0, sigma2 = NA, aicc = NA, se.phi = se.phi,
             se.theta = 0)
    a1 = lapply(1:nrow(Xc), function(i){.innovation.update(Xc[i,], a)})
    a = list(phi = phi, theta = 0,
             sigma2 = do.call(sum, lapply(1:nrow(Xc), function(i){a1[[i]]$sigma2}))/nrow(Xc),
             aicc = do.call(sum, lapply(1:nrow(Xc), function(i){a1[[i]]$aicc}))/nrow(Xc),
             se.phi = se.phi,
             se.theta = 0, p = p)
    return(a)
  }, mc.cores=n_core)
  aicc_threshold = Inf
  for (i in 1:length(all_models)){
    if (all_models[[i]]$aicc < aicc_threshold){
      a = all_models[[i]]
      aicc_threshold = all_models[[i]]$aicc
    }
  }
  return(a)
}



gcica_bss_dwst_single = function(group, M = nrow(group[[1]]), 
                                 W1 = diag(M), tol = 1e-03, 
                                 maxit = 100, nmaxit = 1, maxnmodels = 10, 
                                 prewhite = T, max_stuck = 6,
                                 max_inf_stuck = 3,
                                 iter_print = 10) {

  #################
  # Preprocessing #
  #################
  # M: Number of mixutures (number of sources)
  # N: number of time points
  Xc_rows = nrow(group[[1]])
  if (M > Xc_rows) {
    stop("Number of sources must be less or equal than number of variables")
  }
  N = ncol(group[[1]])

  # Remove subject with different number of columns
  group = group[sapply(group, function(x) ncol(x) == N)]
  
  # Number of subjects in the group
  num_subject = length(group)
  
  # Prewhite
  if (prewhite){
    K = list()
    for (i in 1:num_subject) {
      group[[i]] = t(scale(t(group[[i]]), center=TRUE, scale=FALSE))
      svdcovmat = svd(group[[i]]/sqrt(N))
      K[[i]] = t(svdcovmat$u %*% diag(1/svdcovmat$d))[1:M, ]
      group[[i]] = K[[i]] %*% group[[i]]
      }
  }

  wlik = -Inf

  freqlength = floor(N/2 - 1)
  freq = 1:freqlength * 2 * pi/N

  # g = matrix(0, M, freqlength)

  # Initial source
  WXc = list()
  for (i in 1:num_subject) {
    WXc[[i]] = W1 %*% group[[i]]
  }

  # Discrete Fourier Transformation
  X_dftall = list()
  for (i in 1:num_subject) {
    X_dftall[[i]] = t(mvfft(t(group[[i]])))/sqrt(2 * pi * N)
  }

  # Estimate time series order p and parameters phi
  sourcetsik = lapply(1:M, function(m){
    tmp = matrix(0, nrow = num_subject, ncol = N)
    l = 1
    for (i in 1:num_subject){
      tmp[l, ] = WXc[[i]][m, ]
      l = l + 1
    }
    return(tmp)
  })
  
  g = do.call("rbind", mclapply(1:M, function(m){
    fit = group_yw(sourcetsik[[m]], maxnmodels = maxnmodels)
    if (fit$p == 0){
      return(fit$sigma2/(2 * pi) * rep(1, freqlength))
    }
    else {
      return((fit$sigma2/(2 * pi))/(abs(1 - matrix(fit$phi, 1, fit$p) %*%
                                            exp(-(0+1i) * matrix(1:fit$p, fit$p, 1) %*% freq))^2))
    }
  }, mc.cores=n_core))


  #rm(list = c("WXc"))
  lim = 1
  iter = 0
  NInv = 0
  index1 = as.double(gl(M, M))
  index2 = as.double(gl(M, 1, M^2))
  indx = 2:(freqlength + 1)

  tmp = lapply(1:num_subject, function(j){
    Re(X_dftall[[j]][index1, indx] *
         Conj(X_dftall[[j]][index2, indx]))
  })
  
  stuck = 0
  inf_stuck = 0
  reinit = 0

  while (lim > tol & iter < maxit & NInv < nmaxit) {
    iter = iter + 1
    taucount = 1
    err = 1
    orthoerror = 1
    W2 = W1
    tau = 0.5
    eigenval = rep(0, M)

    while (taucount < 60 & err > 1e-05 & orthoerror > 1e-05) {
      #print(taucount)
      W2 = do.call("rbind", mclapply(1:M, function(m){
        Gam = 0
        if (m > 1) {
          Gam = do.call("sum", lapply(1:(m - 1), function(k){
            matrix(W2[k, ], M, 1) %*% matrix(W2[k, ], 1, M)
          }))
        }
        tmpmat = lapply(1:num_subject, function(j){
          t(matrix(tmp[[j]] %*% matrix(1/g[m, ], freqlength, 1), M, M))})
        
        tmpmat = Reduce('+', tmpmat)
        tmpV = tmpmat + tau * Gam
        eigenv = eigen(tmpV)
        eigenval[m] = eigenv$values[M]
        return(eigenv$vectors[,M])
      }, mc.cores=n_core))
      
      orthoerror = sum(sum((W2 %*% t(W2) - diag(rep(1,M)))^2))
      err = amari_distance(rerow(W1), rerow(W2))
      taucount = taucount + 1
      tau = 2 * tau
    }
    #print(W2)
    wlik2 = -1 * sum(eigenval) - 1 * sum(log(g)) + N * log(abs(det(W2)))
    if (wlik < wlik2) {
      Wbest = W2
      wlik = wlik2
      if (iter %% iter_print == 0){print(paste("Group Color ICA - Iteration", iter,
                  ": update successful."))}
    }
    # else if (wlik2 == -Inf){
    #   W2 = matrix(rnorm(M * M), M, M)
    #   print(paste("Group Color ICA - Iteration", iter,
    #               ": current Whittle likelihood (", wlik2, ") is negative infinity."))
    # }
    else if (wlik2 == -Inf){
      if (iter %% iter_print == 0){print(paste("Group Color ICA - Iteration", iter,
                  ": current Whittle likelihood  is -Inf."))}
      if (inf_stuck >= max_inf_stuck){
        if (iter %% iter_print == 0){print(paste("Whittle likelihood stuck at -Inf for", max_inf_stuck, "times. Reinitialize W2."))}
        W2 = matrix(rnorm(M * M), M, M)
        reinit = reinit + 1
        inf_stuck = 0
        maxit = maxit + inf_stuck
      }
      else {inf_stuck = inf_stuck + 1}
    }
    else {
      if (iter %% iter_print == 0){print(paste("Group Color ICA - Iteration", iter,
                     ": current Whittle likelihood (", wlik2, ") is smaller than the best one (", wlik, ")."))}
      stuck = stuck + 1
      if (stuck >= max_stuck){
        if (iter %% 10 == 0){print(paste("Whittle likelihood stuck for", max_stuck, "times. Reinitialize W2."))}
        W2 = matrix(rnorm(M * M), M, M)
        reinit = reinit + 1
        stuck = 0
        maxit = maxit + stuck
      }
    }
    lim = err
    if (iter %% iter_print == 0){print(paste("Group Color ICA - Iteration", iter, ": error is equal to ",
                lim))}
    W1 = W2
    if ((iter == maxit & NInv < nmaxit)) {
      if (iter %% iter_print == 0){print("Group Color ICA: iteration reaches to maximum. Start new iteration.")}
      W2 = matrix(rnorm(M * M), M, M)
      reinit = reinit + 1
      qrdec = qr(W2)
      W2 = qr.Q(qrdec)
      # iter = 0
      NInv = NInv + 1
    }
    
    WXc = lapply(1:num_subject, function(i){
      W2 %*% group[[i]]
    })
    
    sourcetsik = lapply(1:M, function(m){
      tmp = matrix(0, nrow = num_subject, ncol = N)
      l = 1
      for (i in 1:num_subject){
        tmp[l, ] = WXc[[i]][m, ]
        l = l + 1
      }
      return(tmp)
    })
    
    g = do.call("rbind", mclapply(1:M, function(m){
      fit = group_yw(sourcetsik[[m]])
      if (fit$p == 0){
        return(fit$sigma2/(2 * pi) * rep(1, freqlength))
      }
      else {
        return((fit$sigma2/(2 * pi))/(abs(1 - matrix(fit$phi, 1, fit$p) %*%
                                            exp(-(0+1i) * matrix(1:fit$p, fit$p, 1) %*% freq))^2))
      }
    }, mc.cores=n_core))

    if (NInv == nmaxit) {
      print("Group Color ICA: no convergence")
      }
    }
    if (wlik > wlik2) {
      W2 = Wbest
      wlik2 = wlik
    }


  #wt = lapply(1:num_subject, function(i){Wbest %*% K[[i]]})
  result = new.env()
  result$W = Wbest
  result$wlik = wlik
  #result$K = K
  # result$A = lapply(1:num_subject, function(i){
  #   tryCatch({
  #     t(wt[[i]]) %*% solve(wt[[i]] %*% t(wt[[i]]))
  #   }, error = function(e){NULL})})
  result$S = lapply(1:num_subject, function(i){Wbest %*% group[[i]]})
  result$iter = iter
  result$NInv = NInv
  result$den = g
  result$taucount = taucount
  result$amari = lim
  result$reinit = reinit
  #result$group = group
  result = as.list(result)
  return(result)
}