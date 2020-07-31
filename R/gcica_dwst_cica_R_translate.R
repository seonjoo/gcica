#' Implement different W same T group-colorICA algorithm.
#'
#' @param Xc      Hierachical list containing       [list]
#' data matrices. G is the number of groups and NG is the number of subjects in each group.
#' Each component of the list is a list with data matrices with dimension M by N,
#' where M is the number of mixtures and N is the number of time points.
#' Each data matrix has been pre-whitenned and centered.
#'
#' @return W: a list of M by M unmixing matries
#' @export
#'
#' @examples


gcica_bss_dwst = function (Xc, M = nrow(Xc[[1]][[1]]), W1 = diag(M),
                           tol = 1e-04, maxit = 20, nmaxit = 1,
                           maxnmodels = 100, num_cores = 2) {
  #################
  # Load packages #
  #################
  library(doParallel)
  library(pracma)
  registerDoParallel(cores = num_cores)

  #################
  # Preprocessing #
  #################
  # M: Number of mixutures (number of sources)
  # N: number of time points
  p = nrow(Xc[[1]][[1]])
  if (M > p) {
    stop("Number of sources must be less or equal than number \n  of variables")
  }
  N = ncol(Xc[[1]][[1]])

  # Number of groups
  num_group = length(Xc)

  # Number of subjects in each group
  num_group_subject = lapply(1:num_group, function(i){
    return(length(Xc[[i]]))
  })

  result = foreach(i = 1:num_group) %dopar% {
    # Prewhite
    if (prewhite == 1){
    svdcovmat = Xc
    K = Xc
    for (j in 1:num_group_subject[i]) {
      Xc[[i]][[j]] = t(scale(t(Xc[[i]][[j]]), center=TRUE, scale=FALSE))
      svdcovmat[[i]][[j]] = svd(Xc[[i]][[j]]/sqrt(N))
      K[[i]][[j]] = t(svdcovmat[[i]][[j]]$u %*% diag(1/svdcovmat[[i]][[j]]$d))
      K[[i]][[j]] = K[[i]][[j]][1:M, ]
      Xc[[i]][[j]] = K[[i]][[j]] %*% Xc[[i]][[j]]
    }


    wlik = -Inf

    freqlength = floor(N/2 - 1)
    freq = 1:freqlength * 2 * pi/N


    g = matrix(0, M, freqlength)


    # Initial source
    for (j in 1:num_group_subject[i]) {
      WXc[[i]][[j]] = W1 %*% Xc[[i]][[j]]
    }

    # Discrete Fourier Transformation
    X_dftall = Xc

    for (j in 1:num_group_subject[i]) {
      X_dftall[[i]][[j]] = t(mvfft(t(Xc[[i]][[j]])))/sqrt(2 * pi * N)
    }

    # Estimate time series order p and parameters phi
    # Currently don't know how to dopar
    sourcetsik = rep(NA, M)
    for (m in 1:M) {
      sourcetsik[[m]] = matrix(0, nrow = sum(num_subj_group), ncol = N)
      l = 1
      for (j in 1:num_group_subject[i]){
        l = l + 1
        sourcetsik[[m]][l, ] = WXc[[i]][[j]][m, ]
      }
    }

    for (m in 1:M) {
      fit = ar.yw(sourcetsik[[m]], order.max = maxnmodels)
      if (fit$order == 0){
        g[m, ] = fit$var.pred/(2 * pi) * rep(1, freqlength)
      }
      else {
        g[m, ] = (fit$var.pred/(2 * pi))/(abs(1 - matrix(fit$ar, 1, fit$order) %*%
        exp(-(0+1i) * matrix(1:fit$order, fit$order, 1) %*% freq))^2)
      }
    }

    rm(list = c("WXc"))
    lim = 1
    iter = 0
    NInv = 0
    index1 = as.double(gl(M, M))
    index2 = as.double(gl(M, 1, M^2))
    indx = 2:(freqlength + 1)

    tmp = lapply(1:num_group_subject[i], function(j){
      Re(X_dftall[[i]][[j]][index1, indx] *
           Conj(X_dftall[[i]][[j]][index2, indx]))
    })


    while (lim > tol & iter < maxit & NInv < nmaxit) {
      iter = iter + 1
      taucount = 1
      err = 1
      orthoerror = 1
      W2 = W1
      tau = 0.5
      eigenval = rep(0, M)

      while (taucount < 60 & err > 1e-05 & orthoerror > 1e-05) {
        for (m in 1:M) {
          Gam = 0
          if (m > 1) {
            for (k in 1:(m - 1)) {
              nu = matrix(W2[k, ], M, 1) %*% matrix(W2[k, ], 1, M)
              Gam = Gam + nu
            }
          }

          tmpmat = lapply(1:num_group_subject[i], function(j){
            t(matrix(tmp[[i]][[j]] %*% matrix(1/g[m, ],
                                              freqlength, 1), M, M))
          })

          tmpmat = do.call(sum, tmpmat)
          tmpV = tmpmat + tau * Gam
          eigenv = eigen(tmpV)
          eigenval[m] = eigenv$values[M]
          W2[m, ] = eigenv$vectors[, M]
        }
        orthoerror = sum(sum((W2 %*% t(W2) - diag(rep(1,M)))^2))
        err = amari_distance(rerow(W1), rerow(W2))
        taucount = taucount + 1
        tau = 2 * tau
      }

      wlik2 = -1 * sum(eigenval) - 1 * sum(log(g)) + N * log(abs(det(W2)))
      if (wlik < wlik2) {
        Wtmp = W1
        wlik = wlik2
      }
      else print(paste("Color ICA - Iteration ", iter,
                       ": current Whittle likelihood(", wlik2, ")
                       is smaller than previous one (",
                       wlik, ")."))
      lim = err
      print(paste("Color ICA - Iteration ", iter, ": error is equal to ",
                  lim, sep = ""))
      W1 = W2
      if ((iter == maxit & NInv < nmaxit)) {
        print("Color ICA: iteration reaches to maximum. Start new iteration.")
        W2 = matrix(rnorm(M * M), M, M)
        qrdec = qr(W2)
        W2 = qr.Q(qrdec)
        iter = 0
        NInv = NInv + 1
      }
      WXc = W2 %*% Xc
      for (j in 1:M) {
        fit = ar.yw(WXc[j, ], order.max = maxnmodels)
        if (fit$order == 0)
          g[j, ] = fit$var.pred/(2 * pi) * rep(1, freqlength)
        else g[j, ] = (fit$var.pred/(2 * pi))/(abs(1 -
          matrix(fit$ar, 1, fit$order) %*% exp(-(0+1i) *
          matrix(1:fit$order, fit$order, 1) %*% freq))^2)
      }
      if (NInv == nmaxit) {
        print("Color ICA: no convergence")
      }
    }
    if (wlik > wlik2) {
      W2 = Wtmp
      wlik2 = wlik
    }
    }

    wt = lapply(1:num_group_subject[i], function(j){
      W2 %*% K[[i]][[j]]})
    result = new.env()
    result$W = W2
    result$K = K[[i]]
    result$A = lapply(1:num_group_subject[i], function(j){
      t(wt[[j]]) %*% solve(wt[[j]] %*% t(wt[[j]]))})
    result$S = lapply(1:num_group_subject[i], function(j){
      wt[[j]] %*% Xc[[i]][[j]]})
    result$X = Xc[[i]]
    result$iter = iter
    result$NInv = NInv
    result$den = g
    result = as.list(result)
    return(result)
}
return(result)
  }
