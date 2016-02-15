###########################################################

#------------------fit stats throw error so taken out----------------------
### compute a bunch of fit-related stuff for Bayesian p-values

for(i in 1:nplots){
  y.fit[i,1:8] ~ dmulti(muc[i,1:8,1], ncap[i,1])
  y.fit[i,9:16] ~ dmulti(muc[i,1:8,2], ncap[i,2])
  for(k in 1:2){
    ncap.fit[i,k] ~ dbin(pcap[i,k], N[i,k])
  }
  for(t in 1:8){
    e1[i,t] <- muc[i,t,1] * ncap[i,1]
    resid1[i,t] <- pow(pow(y[i,t],0.5)-pow(e1[i,t],0.5),2)
    resid1.fit[i,t] <- pow(pow(y.fit[i,t],0.5) - pow(e1[i,t],0.5), 2)
    
    e1[i,t+8] <- muc[i,t,2] * ncap[i,2]
    resid1[i,t+8] <- pow(pow(y[i,t+8],0.5) - pow(e1[i,t+8],0.5), 2)
    resid1.fit[i,t+8] <- pow(pow(y.fit[i,t+8],0.5) - pow(e1[i,t+8],0.5), 2)
  }
  
  e2[i,1] <- N[i,1]*lambda[i]
  e2[i,2] <- N[i,2]*gamma[i]
  
  for(k in 1:2){
    resid2[i,k] <- pow(pow(ncap[i,k],0.5) - pow(e2[i,k],0.5),2)
    resid2.fit[i,k] <- pow(pow(ncap.fit[i,k],0.5) - pow(e2[i,k],0.5),2)
  }
}

ft1.data <- sum(resid1[,])
ft1.post <- sum(resid1.fit[,])

ft2.data <- sum(resid2[,])
ft2.post <- sum(resid2.fit[,])