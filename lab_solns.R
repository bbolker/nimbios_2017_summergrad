library(deSolve)

##' discrete-time SIR with observation+process error
SIRdsim_obsproc <- function(t,y,params) {
   g <- with(as.list(c(y,params)), {
      p <- 1-exp(-beta*I/N)
      g <- 1-exp(-gamma)
      newinf <- rbinom(1,size=S,prob=p)
      recover <- rbinom(1,size=I,prob=g)
      c(S=S-newinf,
        I=I+newinf-recover,
        R=R+recover,
        incidence=rbinom(1,size=newinf,prob=rptprob))
  })
   list(g)
}

##' discrete-time deterministic SIR
SIRdsim_determ <- function(t,y,params) {
   g <- with(as.list(c(y,params)), {
      p <- 1-exp(-beta*I/N)
      gg <- 1-exp(-gamma)
      newinf <- p*S
      recover <- gg*I
      ## cat(S,I,newinf,recover,"\n")
      c(S=S-newinf,
        I=I+newinf-recover,
        R=R+recover,
        incidence=newinf)
  })
   list(g)
}

##' generate 'data' for testing

set.seed(101)
m1_obsproc <- ode(t=1:20,
          y=c(S=99,I=1,R=0,incidence=0),
          func=SIRdsim_obsproc,
          parms=c(beta=4.0,gamma=1.0,N=100,rptprob=0.8),
          method="iteration")
saveRDS(m1_obsproc,file="data/sim.rds")

##' trajectory matching
determfun <- function(parms,nobs,N) {
    ## cat(parms,"\n")
    res <- ode(t=1:nobs,
               y=c(S=99,I=1,R=0,incidence=0),
               func=SIRdsim_determ,
               parms=c(parms,N=N),
               method="iteration")
    inc <- res[,"incidence"]
    lines(inc)
    return(inc)
}
determfun(c(beta=4,gamma=1),N=100,nobs=10)
library(bbmle)

plot(m1_obsproc[1:7,"incidence"])
tfit <- mle2(minuslogl=incidence~dpois(determfun(c(beta=beta,gamma=gamma),
                                         nobs=7,N=100)),
     data=list(incidence=m1_obsproc[1:7,"incidence"]),
     trace=TRUE,
     start=list(beta=4,gamma=1),
     method="L-BFGS-B",
     lower=c(0.1,0.1))

lines(determfun(coef(tfit),nobs=7,N=100))
