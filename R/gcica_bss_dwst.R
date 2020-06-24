
#' Implement different W same T group-colorICA algorithm.
#'
#' @param Xc      G by NG cell-array contains observed         [cell-array]
#' data matrices. G is the number of groups and NG is the number of subjects in each group.
#' Each component of the cell-array is a matrix with dimension M by T
#' where M is the number of mixtures and T is the number of time points.
#' Each data matrix has been pre-whitenned and centered.
#'
#' @return W: M by M unmixitng matrix
#' @export
#'
#' @examples
#' x=rnorm(10)
gcica_bss_dwst<-function(Xc){
  print('hello')
  return(W)
}
