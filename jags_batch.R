library(R2jags)
L <- load("data/sim.rda")
incid <- m1_obsproc[,"incidence"]
incid <- incid[incid>0]
jags.inits <- function(){
    list(newinf=incid,beta=4,gamma=1)
}

j1_long <- jags.parallel(model.file="sir.jags",
                data=list(incid=incid,
                          N=100,
                          nobs=length(incid)),
                n.chains=3,
                inits=jags.inits,
                ## n.iter=1000,
                n.iter=100000,
                DIC=FALSE,
                parameters.to.save=c("beta","gamma","p.report"))

saveRDS(j1_long,file="data/j1_long.rds")
