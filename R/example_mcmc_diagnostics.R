library(rjags)
library(coda)
library(mcmc)


data(mcmc)
trace <- mcmc.epi3$trace
head(trace)

plotTraj(data = epi3)
