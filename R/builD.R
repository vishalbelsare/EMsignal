#' Put together inverse covariance matrix of overdifferenced signals.
#'     notated as $D$ in (Livsey and McElroy, 2018+ - EMsigex)
#'
#' @param mdl sigex mdl object
#'
#' @return List of inverse covariance matrices
#' @export
#'

builD = function(mdl){

  TT = length(mdl[[4]][[1]]) # total number observations
  J = length(mdl[[2]]) # total number of components in model
  d = unlist(lapply(mdl[[3]], length)) # diff order of each component
  d.full = sum(d) - 1 # remove irregular diff order of 1

  # put together overdifferencing operator coef vectors
  # NOTE - Jth component is irregular and diff.over[[J]] = diff.full
  diff.over = list() # over differenced component operators
  for(j in 1:J){
    diff.over[[j]] = suppressWarnings(sigex.delta(mdl = mdl, omits = j) )
    diff.over[[j]] = round(diff.over[[j]])
  }

  # toeplitz matrix for each over-differenced signal
  Gam = list()
  for(j in 1:J){
    Gam[[j]] = toeplitz(ARMAacvf(ma = diff.over[[j]], lag.max=(TT-d.full)))
  }

  invGam = lapply(Gam, solve)

  return(invGam)
}
