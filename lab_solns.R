SIRdsim_obsproc <- function(t,y,params) {
   g <- with(as.list(c(y,params)), {
      p <- 1-exp(-beta*dt)
      inf <- rbinom(1,size=S,prob=p)
      recover <- rbinom(1,size=I,prob=gamma*dt)
      c(S=S-inf,
        I=I+inf-recover,
        R=R+recover,
        incidence=rbinom(1,size=inf,prob=rptprob))
  })
   list(g)
}

SIRdsim_determ <- function(t,y,params) {
   g <- with(as.list(c(y,params)), {
      p <- 1-exp(-beta*dt)
      inf <- p*S
      recover <- gamma*dt*I
      c(S=S-inf,
        I=I+inf-recover,
        R=R+recover)
  })
   list(g)
}



