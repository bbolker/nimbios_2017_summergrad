library(deSolve)

##' discrete-time stochastic SIR with observation+process error
SIRdsim_obsproc <- function(t,y,parms) {
   g <- with(as.list(c(y,parms)), {
      prob.inf  <- 1-exp(-beta*I/N)
      prob.rec <- 1-exp(-gamma)
      newinf <- rbinom(1,size=S,prob=prob.inf)
      recover <- rbinom(1,size=I,prob=prob.rec)
      ## cat(S,I,newinf,recover,rptprob,"\n")
      c(S=S-newinf,
        I=I+newinf-recover,
        R=R+recover,
        incidence=rbinom(1,size=newinf,prob=rptprob),
        recover=rbinom(1,size=recover,prob=rptprob))
  })
   list(g)
}

##' discrete-time deterministic SIR
SIRdsim_determ <- function(t,y,parms) {
    ## print(c(y,parms))
    g <- with(as.list(c(y,parms)), {
        prob.inf  <- 1-exp(-beta*I/N)
        prob.rec <- 1-exp(-gamma)
        newinf <- prob.inf*S
        recover <- prob.rec*I
        ## cat(S,I,newinf,recover,"\n")
        c(S=S-newinf,
          I=I+newinf-recover,
          R=R+recover,
          incidence=newinf,
          recover=recover)
    })
   list(g)
}

##' generate 'data' for testing
set.seed(101)
m1_obsproc <- ode(t=1:20,
                  y=c(S=99,I=1,R=0,incidence=1,recover=0),
                  func=SIRdsim_obsproc,
                  parms=c(beta=2.0,gamma=1.0,N=100,rptprob=0.8),
                  method="iteration")
m1_proc <- ode(t=1:20,
               y=c(S=99,I=1,R=0,incidence=1,recover=0),
               func=SIRdsim1,
               parms=c(beta=2.0,gamma=1.0,N=100),
               method="iteration")
m1_obs <- ode(t=1:20,
               y=c(S=99,I=1,R=0,incidence=1,recover=0),
               func=SIRdsim_determ,
               parms=c(beta=2.0,gamma=1.0,N=100),
               method="iteration")
save(list=ls(pattern="m1_.*"),file="data/sim.rda")

##' trajectory matching
##'
##' deterministic function, returns incidence
determfun <- function(parms,nobs=20,N=100) {
    ## cat(parms,"\n")
    res <- ode(t=1:nobs,
               y=c(S=N-1,I=1,R=0,incidence=1,recover=0),
               func=SIRdsim_determ,
               parms=c(parms,N=N),
               method="iteration")
    inc <- res[,"incidence"]
    ## lines(inc)
    return(inc)
}

