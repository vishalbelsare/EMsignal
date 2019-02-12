library(sigex)
library(mvtnorm)
library(EMsigex)

# ---- Simulate Data ----------------------------------------------------------
N = 6
T <- TT <- 300
t = 1:T
Phi=diag(N)
Sig=diag(N); Sig[1,2] <- Sig[2,1] <- .75

s1 = gen_trendComp(T, Phi, Sig)
s2 = gen_seasComp(T, Phi, diag(N))
s0 = rmvnorm(n = T, mean = rep(0,N), sigma = diag(N))

data = s1+s2+s0
data = demean(data)

plot(ts(data), main="simulated series")

# ---- Modeling ---------------------------------------------------------------
transform = "none"
x <- t(data)
N <- dim(x)[1]
T <- TT <- dim(x)[2]

# ---- Model ------------------------------------------------------------------
def <- c(0,1,0,1)
mdl = NULL
mdl = sigex.add(mdl, seq(1,N), 'wn', c(1,-1), def)    # Trend
mdl = sigex.add(mdl, seq(1,N), 'wn', rep(1, 12), def) # Seasonal
mdl <- sigex.add(mdl,seq(1,N),"wn",1,def)             # Irregular
mdl <- sigex.meaninit(mdl,data,0)

# Set default parameters
par.default <- sigex.default(mdl,data)[[1]]
flag.default <- sigex.default(mdl,data)[[2]]
psi.default <- sigex.par2psi(par.default,flag.default,mdl)

# Set param to TRUE values
param = par.default
