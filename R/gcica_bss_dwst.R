
#' Implement different W same T group-colorICA algorithm.
#'
#' @param Xc      Hierachical list containing       [list]
#' data matrices. G is the number of groups and NG is the number of subjects in each group.
#' Each component of the list is a list with data matrices with dimension M by N,
#' where M is the number of mixtures and N is the number of time points.
#' Each data matrix has been pre-whitenned and centered.
#'
#' @return W: M by M unmixitng matrix
#' @export
#'
#' @examples
#' x=rnorm(10)
gcica_bss_dwst = function(Xc, scovmatdnmntropt = 'sz', prewhite = 1, M = nrow(Xc[[1, 1]]),
                          Win = diag(M), arma = 0, tol = 1e-10, maxit = 200, num_initials = 10,
                          num_cores = 2){
  #################
  # Load packages #
  #################
  library(doParallel)
  library(pracma)
  registerDoParallel(cores = num_cores)
  #################
  # Preprocessing #
  #################
  # M: Number of  mixutures (number of sources)
  # N: number of time points
  M = nrow(Xc[[1]][[1]])
  N = ncol(Xc[[1]][[1]])

  # Number of groups
  num_group = length(Xc)

  # Number of subjects in each group
  num_group_subject = lapply(1:num_group, function(i){
    return(length(Xc[[i]]))
  })

  #  Denominator used to calculate sample covariance matrix
  #       'sz': using sample size N [default]
  #       'szm1': sample size - 1 (N - 1)
  if (scovmatdnmntropt == 'sz'){
    scovmatdnmntr = N
  } else{
    scovmatdnmntr = N - 1
  }


  # Prewhite
  if (prewhite == 1){
    Xc = foreach(i = 1:num_group) %dopar% {
      Xc_i = foreach(j = 1:num_group_subject[i]) %dopar% {
        Xc[[i]][[j]] = t(scale(t(Xc[[i]][[j]]), center=TRUE, scale=FALSE))
        svd.fit = svd(Xc[[i]][[j]]) #Xc = UDV’
        Xprew = Xc[[i]][[j]] %*% svd.fit$v ## n x p matrix
        # Xc[[i]][[j]] = Xprew
        return(Xprew)
      }
      return(Xc_i)
      }
    }


  # Discrete Fourier Transformation
  X_dftall =
  # matrix(data = NA, nrow = nrow(Xc), ncol = ncol(Xc))
  foreach (i = 1:num_group) %dopar% {
    X_dft_i = foreach (j = 1:num_group_subject[i]) %dopar% {
      X_dftall[[i, j]] = t(mvfft(t(Xc[[i, j]])))/sqrt(2 * pi * N)
      return(X_dftall[[i, j]])
    }
    return(X_dft_i)
  }

  # Frequencies
  freq = t(0:N-1)/N*2*pi

  ###############################################
  # Initial W and source time series parameters #
  ###############################################

  # Initial W
  W1 = Win
  # Initial source
  WXc = Xc
  WXc = foreach (i = 1:num_group) %dopar% {
    wX_i = foreach (j = 1:num_group_subject[i]) %dopar% {
      WXc[[i]][[j]] = W1[i] %*% Xc[[i]][[j]]
      return(WXc[[i]][[j]])
    }
    return(wX_i)
  }

  # Estimate time series order p and parameters phi
  # sourcetsik = matrix(data = NA, nrow = 1, ncol = M)
  sourcetsik = foreach (m = 1:M) %dopar% {
    sourcetsik[m] = matrix(0, nrow =sum(num_subj_group), ncol = N)
    l = 1
    for (i in 1:num_group) %dopar% {
      for (j in 1:num_group_subject[i]) %dopar% {
        l = l + 1
        sourcetsik[m][l, ] = wX[[i, j]][m, ]
      }
    }
  }

  for (m in 1:M){
    if (arma == 1){
      stop('We only consider AR model for now!')
    }
    if (arma == 0){
      fit = ar.yw(sourcetsik[m], order.max = maxnmodels)
      if (fit$order == 0){
        g[m, ] = fit$var.pred/(2 * pi) * rep(1, N)
      }
      else {
        g[m, ] = (fit$var.pred/(2 * pi))/
          (abs(1 - matrix(fit$ar, 1, fit$order) %*%
                 exp(-(0+1i) * matrix(1:fit$order, fit$order, 1) %*% freq))^2)
      }
    }
    else {
      stop('Wrong arma value!')
    }
    p[m] = fit$order
    sigmasq[m] = fit$var.pred
    phi[m] = iphi;
  }


  ################################
  # Iterative updating algorithm #
  ################################

  # Distance between W1 and W2 (previous and current iteration)
  wddist = 1
  # Loop control
  iter = 0
  # How many initial W tried
  num_initials = 1
  # How many singular Hessian matrix found
  num_sg = 0
  lambda = repmat(matrix(1, nrow = M*(M+1)/2, ncol = 1), 1, num_group)
  index1 =  reshape(reshape(t(repmat(1:M,1,M),M,M)), 1, M^2);
  index2 = repmat(1:M,1,M)
  tempmx11 = reshape(index2,M,M)
  tempmx22 = reshape(index1,M,M)
  # row indices for Lower triangular elements
  tempmx1 = tempmx11[index1<=index2]
  # column indices for lower triangular elements
  tempmx2 = tempmx22[index1<=index2]

  while (wddist > tol & iter < maxit + 1 & NInv < maxit + 1){
    #################################################################
    # Newton-Raphson method with Lagrange multiplier for updating W #
    #################################################################
    iter = iter + 1
    for (i in 1:num_group){
      for (j in 1:num_group_subject[i]){
        X_dft[[i, j]] = X_dftall[[i, j]]
      }
      W1_inv[i, ] = inv(W1[i, ])
    }

    for (i in 1:num_group){
      # Objective function
      # F = -L + CC
      # Goal: minimize F
      # -L = L1 + L2
      # L1 = \frac{1}{2} sum_j sum_k{} in Eqn. (5) in JASA 2011
      # L2 = -T \log|det(W)|
      # CC = \lambda'*C
      # Let AA = WW^T - I. Then C(1) is AA(1,1), C(2) is AA(2,1), C(3)
      # is AA(3,1), ..., C(M) is AA(M,1), C(M+1) is AA(2,2), ... C
      # contains the lower diagonal part of AA and goes from up to down
      # and left to right.
      #
      # First derivative for score function
      # Take derivative w.r.t (W_11,
      # W_12,...,W_1M,W_21,...,W_2M,...,W_M1,...W_MM)
      #
      # dF/ d(W_11, W_12, W_13,...,W_MM)
      Score_l1_w = matrix(0, M^2, 1)
      Score_l2_w = matrix(0, M^2, 1)
      for (j in 1:num_group_subject[i]){
        Score_l1_w = Score_l1_w + sum(real(conj(X_wdft{ik}{jk}(index1,:)) .* X_dft{ik}{jk}(index2,:))./ ...
                                      g(index1,:),2)
        Score_l2_w = Score_l2_w - T * reshape( W1_inv{ik}, M^2, 1)
      }
    }
  }

  svdcovmat = svd(Xc/sqrt(N))
  K = t(svdcovmat$u %*% diag(1/svdcovmat$d))
  K = K[1:M, ]
  Xc = K %*% Xc
  W1 = Win
  wlik = -Inf
  rm(list = c("Win"))
    freqlength = floor(N/2 - 1)
    freq = 1:freqlength * 2 * pi/N
    g = matrix(0, M, freqlength)
    X.dft = t(mvfft(t(Xc)))/sqrt(2 * pi * N)
    WXc = W1 %*% Xc
    for (j in 1:M) {
      fit = ar.yw(WXc[j, ], order.max = maxnmodels)
      if (fit$order == 0)
        g[j, ] = fit$var.pred/(2 * pi) * rep(1, freqlength)
      else g[j, ] = (fit$var.pred/(2 * pi))/(abs(1 - matrix(fit$ar,
                                                            1, fit$order) %*% exp(-(0+1i) * matrix(1:fit$order,
                                                                                                   fit$order, 1) %*% freq))^2)
    }
    rm(list = c("WXc"))
    lim = 1
    iter = 0
    NInv = 0
    index1 = as.double(gl(M, M))
    index2 = as.double(gl(M, 1, M^2))
    indx = 2:(freqlength + 1)
    tmp = Re(X.dft[index1, indx] * Conj(X.dft[index2, indx]))
    while (lim > tol & iter < maxit & NInv < nmaxit) {
      iter = iter + 1
      taucount = 1
      err = 1
      orthoerror = 1
      W2 = W1
      tau = 0.5
      eigenval = rep(0, M)
      while (taucount < 60 & err > 1e-05 & orthoerror >
             1e-05) {
        for (j in 1:M) {
          Gam = 0
          if (j > 1) {
            for (k in 1:(j - 1)) {
              nu = matrix(W2[k, ], M, 1) %*% matrix(W2[k,
                                                       ], 1, M)
              Gam = Gam + nu
            }
          }
          tmpmat = t(matrix(tmp %*% matrix(1/g[j, ],
                                           freqlength, 1), M, M))
          tmpV = tmpmat + tau * Gam
          eigenv = eigen(tmpV)
          eigenval[j] = eigenv$values[M]
          W2[j, ] = eigenv$vectors[, M]
        }
        orthoerror = sum(sum((W2 %*% t(W2) - diag(rep(1,
                                                      M)))^2))
        err = amari_distance(rerow(W1), rerow(W2))
        taucount = taucount + 1
        tau = 2 * tau
      }
      wlik2 = -1 * sum(eigenval) - 1 * sum(log(g)) + N *
        log(abs(det(W2)))
      if (wlik < wlik2) {
        Wtmp = W1
        wlik = wlik2
      }
      else print(paste("Color ICA - Iteration ", iter,
                       ": current Whittle likelihood(", wlik2, ") is smaller than previous one (",
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

  wt = W2 %*% K
  result = new.env()
  result$W = W2
  result$K = K
  result$A = t(wt) %*% solve(wt %*% t(wt))
  result$S = wt %*% Xin
  result$X = Xin
  result$iter = iter
  result$NInv = NInv
  result$den = g
  as.list(result)
}
