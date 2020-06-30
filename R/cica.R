#' To extract source signals from a mixture X. X is an M by T matrix
#' where M index the voxles and T index time points.
#'
#' @param X       M by T matrix. M is the number of observed         [matrix]
#' mixtures (also the number of sources) and T is the
#' number of time series points. X will be centered
#' within the algorithm.
#'
#' @return W: M by M unmixitng matrix
#' @export
#'
#' @examples


cICA = function (Xin, M = dim(Xin)[1], Win = diag(M), tol = 1e-04, maxit = 20,
                 nmaxit = 1, unmixing.estimate = "eigenvector", maxnmodels = 100, batch_size = 32) {
  p = dim(Xin)[1]
  if (M > p) {
    stop("Number of sources must be less or equal than number \n  of variables")
  }
  if (unmixing.estimate != "eigenvector" && unmixing.estimate !=
      "newton") {
    stop("Methods to estimate the unmixing matrix can be \n  'eigenvector' or 'newton' only")
  }
  N = ncol(Xin)
  Xc = t(scale(t(Xin), center = TRUE, scale = FALSE))
  svdcovmat = svd(Xc/sqrt(N))
  K = t(svdcovmat$u %*% diag(1/svdcovmat$d))
  K = K[1:M, ]
  Xc = K %*% Xc
  W1 = Win
  wlik = -Inf
  rm(list = c("Win"))
  if (unmixing.estimate == "eigenvector") {
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
