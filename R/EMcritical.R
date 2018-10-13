#' Calculate updated Sigma values at EM critical point
#'
#' @param j component of mdl
#' @param lMS list with lMS[[1]] = M an [N,N,T,T] array and lMS[[2]] = S[T,N]
#' matrix
#' @param Sig list with current Sigma parameter estimates
#' @param mdl sigex model object
#' @param invGam list of inverse covariance matrices for over-diff signals
#'
#' @return NxN updated Sigma matrix
#' @export
#'

EMcritical = function(j, Sig, lMS, mdl, invGam){

  N = dim(Sig[[j]])[1]
  TT = length(mdl[[4]][[1]])

  d = unlist(lapply(mdl[[3]], length)) - 1 # diff order of each component
  d = sum(d) # full differencing order

  Mj = lMS[[1]][[j]]
  Mj = block2array(Mj, N, TT-d-1)
  S = lMS[[2]]

  # put together overdifferencing operator coef vectors
  # diff.over = sigex.delta(mdl = mdl, omits = j)

  # Mj = M[[j]][,,k,el] - M[[j]][,,k-12,el] -
  #      M[[j]][,,k,el-12] + M[[j]][,,k-12,el-12]

  invGam = builD(mdl = mdl) # THIS IS GETTING CALLED TOO MANY TIMES! FIX LATER

  outSig = matrix(0, N, N)
  for(k in (d+1):(TT-1)){
    for(el in (d+1):(TT-1)){

      # print(sprintf("k = %i   el = %i", k-d, el-d))

      new.term =  invGam[[j]][k-d, el-d] * (Mj[, k-d, , el-d] + (S[[j]][k-d, ] %*% t(S[[j]][el-d, ])))
      outSig = outSig + new.term
    }}
  outSig = outSig / (TT-d)
  return(outSig)
  #return(( outSig + t(outSig) )/2)
}
