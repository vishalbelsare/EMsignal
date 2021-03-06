#' Runs sigex for EM algorithm
#'
#' @param param parameter list in sigex form
#' @param data data matrix (dim x sample_size)
#' @param mdl sigex mdl
#'
#' @return list with M (error matricies) and S (signal estimates)
#' @export
#'

sigexRun = function(param, data, mdl){

  N = dim(data)[2]
  TT = dim(data)[1]

  J = length(mdl[[2]]) # total number of components in model
  d = lapply(mdl[[3]], length) # diff order of each component

  signal = rep(list(matrix(-99, N*TT, N*TT)), J) # initalize storage
  for(j in 1:J){
    signal[[j]] = sigex.signal(data, param, mdl, j)
  }

  extract = rep(list(matrix(-99, TT, N)), J) # initalize storage
  for(j in 1:J){
    extract[[j]] = sigex.extract(data, signal[[j]], mdl, param)
  }

  # put together full differencing operator
  diff.full = sigex.delta(mdl = mdl, omits = NA) # full difference operator
  diff.full = round(diff.full) # this is not needed in updated sigex package
  d.full = length(diff.full)


  # form output into paper notation
  # M = list()
  M.diff = list()
  for(j in 1:J){

    print(j)

    # M[[j]] = block2array(signal[[j]][[2]], N = N, TT = TT)
    # # need to define d the full differencing order
    # M.diff[[j]] = array(NA, c(N, TT-d.full, N, TT-d.full))

    M.diff[[j]] = matrixDiff(m = signal[[j]][[2]], N = N, TT = TT, delta = diff.full)

  }


  S = list()
  for(j in 1:J){
    S[[j]] = extract[[j]][[1]]
    S[[j]] = D %*% S[[j]]
  }

  return(list(M, S))
}
