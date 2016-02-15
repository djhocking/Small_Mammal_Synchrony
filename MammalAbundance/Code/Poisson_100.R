##########--------Dail-Madsen Model with binomial series 100m buffer---------------
# Model
cat("
  model{
    # Abundance Priors
    a.N ~ dnorm(0, 0.01)
    b.elev ~ dnorm(0, 0.01)
    #b.elev2 ~ dnorm(0, 0.01)
    #b.slope ~ dnorm(0, 0.01)
    b.hard100 ~ dnorm(0, 0.01)
    b.mix100 ~ dnorm(0, 0.01)
    #b.soft100 ~ dnorm(0, 0.01)
    b.age100 ~ dnorm(0, 0.01)
    b.stream100 ~ dnorm(0, 0.01)
    #b.wet100 ~ dnorm(0, 0.01)
    b.trap ~ dnorm(0, 0.01)
    b.year2 ~ dnorm(0, 0.01)
    b.year3 ~ dnorm(0, 0.01)
    gam1 ~ dunif(-1, 1)
    gam.a ~ dnorm(0, 0.01)
    
    #omega ~ dunif(0, 1)
    
    # Random priors
    for(i in 1:nsites){
      b.site[i] ~ dnorm(0, tau.site)
      gam0[i] ~ dnorm(0, tau.gam0)
    }
    
    # Hyperpriors for random effects
    sigma.site ~ dunif(0, 5)
    tau.site <- 1/(sigma.site*sigma.site)
    
    sigma.gam0 ~ dunif(0, 5)
    tau.gam0 <- 1/(sigma.gam0*sigma.gam0)
    
    # Detection Priors
    p0 ~ dunif(-5, 5)
    #logitp0 <- log(p0 / (1-p0))
    p.precip ~ dunif(-5, 5)
    p.trap ~ dunif(-5, 5)
    p.temp ~ dunif(-5, 5)
    p.temp2 ~ dunif(-5, 5)
    # p.day ~ dunif(-5, 5)
    
    
    # Abundance
    for(k in 1:1){
      for(i in 1:nplots){    
        # Abundance for closed period 1
        log(lambda[i]) <- a.N + b.site[site[i]] + b.trap*traptype[i] + b.year2*year2[i] + b.year3*year3[i] + b.elev*elev[i] +  b.age100 * age.100[i] + b.stream100 * stream.100[i] + b.mix100 * mix.100[i] + b.hard100 * hard.100[i] # b.wet100 * wetland.100[i] + 
        #z[i] ~ dbern(omega)
        
        #eff.lam[i] <- z[i] * lambda[i] + 0.00001 # hack for JAGS based on Bolker/NCEAS non-linear working group
        
        N[i,k] ~ dpois(lambda[i])
        
        # Open Period Change in Abundance
        gamma[i] <- min( exp(gam.a + gam0[site[i]] + gam1* (N[i,1] - sum(y[i, ,1]))), 1000) # add year2 and year3 effects and time between periods or day of the year. Also consider putting the random effect outside the exponential
        
        #eff.gam[i] <- z[i] * gamma[i] + 0.00001
        N[i,2] ~ dpois(gamma[i])
      }
    }
    
    # Detection - binomial series rather than multinomial formulation
    for(k in 1:2){
      for(i in 1:nplots){
        fac[i,1,k] <- 0
        y[i,1,k] ~ dbin(p[i,1,k], N[i,k] - fac[i,1,k])
        lp[i,1,k] <- p0 + p.precip*precip[i,1,k] + p.trap*traptype[i] + p.temp*temp[i,1,k] + p.temp2*temp[i,1,k]*temp[i,1,k]
        lp.lim[i,1,k] <- min(999, max(-999, lp[i,1,k]))
        p[i,1,k] <- 1/(1 + exp(-lp.lim[i,1,k]))
        
        for(j in 2:8) {
          fac[i,j,k] <- fac[i,j-1,k] + y[i,j-1,k]
          y[i,j,k] ~ dbin(p[i,j,k], N[i,k] - fac[i,j,k])
          lp[i,j,k] <- p0 + p.precip*precip[i,j,k] + p.trap*traptype[i] + p.temp*temp[i,j,k] + p.temp2*temp[i,j,k]*temp[i,j,k]
          lp.lim[i,j,k] <- min(999, max(-999, lp[i,j,k]))
          p[i,j,k] <- 1/(1 + exp(-lp.lim[i,j,k]))
        }
      }
    }
    
    # Fit statistics
    for(k in 1:2) {
      for(i in 1:nplots) {
        fac2[i,1,k] <- 0
        for(j in 1:1) {
          # Compute fit statistic for observed data
          count.hat[i,j,k] <- round((N[i,k] - fac2[i,j,k]) * p[i,j,k]) # predicted
          error.obs[i,j,k] <- y[i,j,k] - count.hat[i,j,k]
          E.obs[i,j,k] <- pow((y[i,j,k] - count.hat[i,j,k]), 2) / (count.hat[i,j,k] + 0.5)
          # Generate replicate data and compute fit stats for them
          count.new[i,j,k] ~ dbin(p[i,j,k], N[i, k] - fac2[i,j,k])
          E.new[i,j,k] <- pow((count.new[i,j,k] - count.hat[i,j,k]),2) / (count.hat[i,j,k] + 0.5)
        }
      for(j in 2:8) {
          fac2[i,j,k] <- fac2[i,j-1,k] + count.hat[i,j-1,k]
          # Compute fit statistic for observed data
          count.hat[i,j,k] <- round((N[i,k] - fac2[i,j,k]) * p[i,j,k]) # predicted
          error.obs[i,j,k] <- y[i,j,k] - count.hat[i,j,k]
          E.obs[i,j,k] <- pow((y[i,j,k] - count.hat[i,j,k]), 2) / (count.hat[i,j,k] + 0.5)
          # Generate replicate data and compute fit stats for them
          count.new[i,j,k] ~ dbin(p[i,j,k], N[i,k] - fac2[i,j,k])
          E.new[i,j,k] <- pow((count.new[i,j,k] - count.hat[i,j,k]),2) / (count.hat[i,j,k] + 0.5)
        }
      }
    }
    
    fit.obs <- sum(E.obs[,,])
    fit.new <- sum(E.new[,,])
    
  }", fill = TRUE, file = 'Model/dm100p.txt')


#------Function to package and run Poisson 100m model-------

dm100p <- function(y, n.burn = 1, n.it = 3000, n.thin = 1, outfile = "Diagnostic_Plots.pdf", outfile2 = "Table.csv") {
  force(n.burn); force(n.it); force(n.thin)
  Name <- y
  nplots <- nrow(y)
  
  y <- array(c(as.matrix(y[,2:9]), as.matrix(y[ ,10:17])), c(length(y$Plot), 8, 2))
  precip.array <- array(c(as.matrix(precip.s[,1:8]), as.matrix(precip.s[ ,9:16])), c(dim(precip.s)[1], 8, 2))
  temp.array <- array(c(as.matrix(Temp.s[,1:8]), as.matrix(Temp.s[ ,9:16])), c(dim(Temp.s)[1], 8, 2))
  
  ymax1 <- apply(y[,1:8, 1],1,sum)
  ymax1 <- ymax1[!is.na(ymax1)]
  ymax1 <- rep(ymax1, length.out = nplots)
  ymax2 <- apply(y[,1:8, 2],1,sum)
  ymax2 <- ymax2[!is.na(ymax2)]
  ymax2 <- rep(ymax2, length.out = nplots)
  ncap <- as.matrix(cbind(ymax1, ymax2))
  
  # Bundle data for JAGS/BUGS
  jdata100 <- list(y=y, nplots=nplots, ncap=ncap, 
                   elev=elev.s,
                   #slope=slope.s,
                   site = as.numeric(Habitat$Site),
                   year2=year2,
                   year3= year3,
                   #Day = as.matrix(Day[ , 2:17]),
                   traptype=as.numeric(Trap[,3]),
                   precip=precip.array,
                   temp = temp.array,
                   #wet.100 = Landscape.s[ , "Wetland_100"]
                   # dist.wet = Landscape.s[, "D_wetland_buff"],
                   age.100 = Landscape.s[ , "Age_100"],
                   hard.100 = Landscape.s[ , "H_100"],
                   mix.100 = Landscape.s[ , "M_100"],
                   stream.100 = Landscape.s[ , "Stream_100"],
                   #day=as.matrix(day.s),
                   nsites = length(unique(Habitat$Site)),
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
         b.site=rnorm(length(unique(Habitat$Site)), 0, 1))
  }
  
  # Set parameters of interest to monitor and save
  params100 <- c("N",
                 "p0",
                 "p.precip",
                 "p.temp",
                 "p.temp2",
                 #"p.day",
                 "p.trap",
                 "a.N",
                 "b.elev",
                 # "b.elev2",
                 "b.mix100",
                 #"b.slope",
                 "b.hard100",
                 "b.age100",
                 "b.stream100",
                 "b.trap",
                 "b.year2",
                 "b.year3",
                 "sigma.gam0",
                 "gam.a",
                 "gam1",
                 #"z",
                 #"omega",
                 "sigma.site", 
                 "count.hat",
                 "fit.obs",
                 "fit.new",
                 "p")
  
  # Run JAGS in parallel for improved speed
  CL <- makeCluster(3) # set number of clusters = to number of desired chains
  clusterExport(cl=CL, list("jdata100", "params100", "inits100", "ymax1", "ymax2", "n.burn", "n.it", "n.thin", "Habitat"), envir = environment()) # make data available to jags in diff cores
  clusterSetRNGStream(cl = CL, iseed = 5312)
  
  out <- clusterEvalQ(CL, {
    library(rjags)
    load.module('glm')
    jm <- jags.model("Model/dm100p.txt", jdata100, inits100, n.adapt = n.burn, n.chains = 1)
    fm <- coda.samples(jm, params100, n.iter = n.it, thin = n.thin)
    return(as.mcmc(fm))
    
  })
  
  out.list <- mcmc.list(out) # group output from each core into one list
  stopCluster(CL)
  
  library(ggmcmc)
  ggmcmc(ggs(out.list[ , c("p0", "p.precip", "p.temp", "p.temp2", "p.trap", "a.N", "b.trap", "b.mix100", "b.hard100", "b.age100", "b.stream100", "b.year2", "b.year3", "b.elev", "gam.a", "sigma.gam0", "gam1")]), file = outfile)
  
  foo <- summary(out.list[ , c("p0", "p.precip", "p.temp", "p.temp2", "p.trap", "a.N", "b.trap", "b.mix100", "b.hard100", "b.age100", "b.stream100", "b.year2", "b.year3", "b.elev", "gam.a", "sigma.gam0", "gam1")])
  
  bar <- data.frame(as.character(row.names(foo[[1]])), foo$statistics[ , "Mean"], foo$quantiles[ , c("2.5%", "97.5%")])
  names(bar) <- c("Variable", "Mean", "2.5%", "97.5%")
  
  write.table(bar, file = outfile2, row.names=F, sep=",")
  
  return(out.list)
}

