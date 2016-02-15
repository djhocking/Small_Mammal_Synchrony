##########--------Dail-Madsen Model with binomial series 100m buffer w/ZIP---------------
# Model
cat("
    model{
    # Abundance Priors
    a.N ~ dnorm(0, 0.01)
    b.elev ~ dnorm(0, 0.01)
    #b.elev2 ~ dnorm(0, 0.01)
    #b.slope ~ dnorm(0, 0.01)
    b.hard100 ~ dnorm(0, 0.01)
    b.soft100 ~ dnorm(0, 0.01)
    b.age100 ~ dnorm(0, 0.01)
    b.stream100 ~ dnorm(0, 0.01)
    b.trap ~ dnorm(0, 0.01)
    b.year2 ~ dnorm(0, 0.01)
    b.year3 ~ dnorm(0, 0.01)
    gam1 ~ dunif(-10, 1)
    #gam.a ~ dnorm(0, 0.01)
    
    # Random priors
    for(i in 1:nsites){
    b.site[i] ~ dnorm(0, tau.site)
    gam0[i] ~ dnorm(0, tau.gam0)
    }
    
    # Hyperpriors for random effects
    sigma.site ~ dunif(0, 3)
    tau.site <- 1/(sigma.site*sigma.site)
    
    sigma.gam0 ~ dunif(0, 3)
    tau.gam0 <- 1/(sigma.gam0*sigma.gam0)
    
    for(i in 1:108){
    for(j in 1:8){
    for(k in 1:2){
    p.delta[i,j,k] ~ dnorm(0, tau.delta)
    }
    }
    }
    
    sigma.delta ~ dunif(0, 5)
    tau.delta <- 1/(sigma.delta*sigma.delta)
    
    
    # Detection Priors
    p0 ~ dunif(-5, 5)
    #logitp0 <- log(p0 / (1-p0))
    p.precip ~ dunif(-5, 5)
    p.trap ~ dunif(-5, 5)
    # p.day ~ dunif(-5, 5)
    
    # Abundance
    for(k in 1:1){
    for(i in 1:nplots){    
    # Abundance for closed period 1
    log(lambda[i]) <- a.N + b.site[site[i]] + b.trap*traptype[i] + b.year2*year2[i] + b.year3*year3[i] + b.elev*elev[i] +  b.age100 * age.100[i] + b.stream100 * stream.100[i] + b.soft100 * soft.100[i] #+ b.hard100 * hard.100[i]
    
    N[i,k] ~ dpois(lambda[i])
    
    # Open Period Change in Abundance
    gamma[i] <- min( exp(gam0[site[i]] + gam1* (N[i,1] - sum(y[i, ,1]))), 500) # gam.a +
    
    N[i,2] ~ dpois(gamma[i])
    }
    }
    
    # Detection - binomial series rather than multinomial formulation
    for(k in 1:2){
    for(i in 1:nplots){
    fac[i,1,k] <- 0
    y[i,1,k] ~ dbin(p[i,1,k], N[i,k] - fac[i,1,k])
    lp[i,1,k] <- p0 + p.precip*precip[i,1] + p.trap*traptype[i] + p.delta[i,1,k]
    lp.lim[i,1,k] <- min(999, max(-999, lp[i,1,k]))
    p[i,1,k] <- 1/(1 + exp(-lp.lim[i,1,k]))
    
    for(j in 2:8) {
    fac[i,j,k] <- fac[i,j-1,k] + y[i,j-1,k]
    y[i,j,k] ~ dbin(p[i,j,k], N[i,k] - fac[i,j,k])
    lp[i,j,k] <- p0 + p.precip*precip[i,j] + p.trap*traptype[i] + p.delta[i,j,k]
    lp.lim[i,j,k] <- min(999, max(-999, lp[i,j,k]))
    p[i,j,k] <- 1/(1 + exp(-lp.lim[i,j,k]))
    }
    }
    }
    }", fill = TRUE, file = 'dm100pod.txt')

#------Function to package and run zip 100m model-------
dm100pod <- function(y, n.burn = 1, n.it = 3000, n.thin = 1) {
  force(n.burn); force(n.it); force(n.thin)
  y <- array(c(as.matrix(y[,2:9]), as.matrix(y[ ,10:17])), c(length(y$Plot), 8, 2))
  nplots <- nrow(y)
  
  ymax1 <- apply(y[,1:8, 1],1,sum)
  ymax2 <- apply(y[,1:8, 2],1,sum)
  ncap <- as.matrix(cbind(ymax1, ymax2))
  
  Precip1 <- Precip
  Precip1[is.na(Precip)] <- 0
  
  # Bundle data for JAGS/BUGS
  jdata100 <- list(y=y, nplots=nplots, ncap=ncap, 
                   elev=elev.s,
                   #slope=slope.s,
                   site=as.factor(Habitat$Site),
                   year2=year2,
                   year3= year3,
                   Day = as.matrix(Day[ , 2:17]),
                   traptype=as.numeric(Trap[,3]),
                   precip=as.matrix(Precip1[,2:9]),
                   # dist.wet = Landscape.s[, "D_wetland_buff"],
                   age.100 = Landscape.s[ , "Age_100"],
                   hard.100 = Landscape.s[ , "H_100"],
                   soft.100 = Landscape.s[ , "S_100"],
                   stream.100 = Landscape.s[ , "Stream_100"],
                   #day=as.matrix(day.s),
                   nsites=nplots/2,
                   nyears=3)
  
  # Set initial values for Gibbs sampler
  inits100 <- function(){
    list(p0=runif(1, 1.1, 2),
         #p.precip=runif(1, 0, 0.1),
         #p.day = runif(1, -.5, 0.1),
         #p.trap=rnorm(1, 0, 0.11),
         #a.N=runif(1,1,1), 
         b.elev=rnorm(1,0,1),
         # b.elev2=rnorm(1,0,1), 
         N=as.matrix(cbind(ymax1, ymax2)), 
         gam1 = rnorm(1, 0.1, 0.1),
         #gam0 = rnorm(54, 0, 0.01),
         #b.trap = runif(1, 0.2, 1),
         b.site=rnorm(54, 0, 1))
  }
  
  # Set parameters of interest to monitor and save
  params100 <- c("N",
                 "p0",
                 "p.precip",
                 #"p.day",
                 "p.trap",
                 "a.N",
                 "b.elev",
                 # "b.elev2",
                 "b.soft100",
                 #"b.slope",
                 #"b.hard100",
                 "b.age100",
                 "b.stream100",
                 "b.trap",
                 "b.year2",
                 "b.year3",
                 #"gam.a",
                 "sigma.gam0",
                 "gam1",
                 "sigma.site",
                 "sigma.delta")
  
  # Run JAGS in parallel for improved speed
  CL <- makeCluster(3) # set number of clusters = to number of desired chains
  clusterExport(cl=CL, list("jdata100", "params100", "inits100", "ymax1", "ymax2", "n.burn", "n.it", "n.thin"), envir = environment()) # make data available to jags in diff cores
  clusterSetRNGStream(cl = CL, iseed = 5312)
  
  out <- clusterEvalQ(CL, {
    library(rjags)
    load.module('glm')
    jm <- jags.model("dm100pod.txt", jdata100, inits100, n.adapt = n.burn, n.chains = 1)
    fm <- coda.samples(jm, params100, n.iter = n.it, thin = n.thin)
    return(as.mcmc(fm))
    
  })
  
  out.list <- mcmc.list(out) # group output from each core into one list
  stopCluster(CL)
  
  #print(plot(out.list[ , c("p0", "p.precip","p.trap","a.N", "b.soft100", "b.hard100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "gam0", "gam1", "sigma.site")]))
  
  return(out.list)
}