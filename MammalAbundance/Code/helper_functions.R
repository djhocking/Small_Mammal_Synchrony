
abund <- function(model = dm100p.NAIN, nplots = nrow(NAIN), output = "Output/Abundance.csv", sep = ",") {
  N.sp <- matrix(NA, nplots, 6)
  for(i in 1:nplots){
    for(k in 1:2){
      foo <- apply(as.matrix(model[, c(paste("N[", i,",", k, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.5, 0.025, 0.975))
      foo <- as.integer(foo)
      for(p in 1:3){
        if(k == 1) {
          N.sp[i,p] <- foo[p]
        }
        if(k == 2){
          N.sp[i,p+3] <- foo[p]
        }
      }
    }
  }
  N.sp <- as.data.frame(N.sp)
  names(N.sp) <- c("N1", "LCI_1", "UCI_1", "N2", "LCI_2", "UCI_2")
  N.sp$Plot <- NAIN$Plot
  write.table(x=N.sp, file=output, row.names=FALSE, sep = ",")
  return(N.sp)
}


# adjust below so it takes the mean for each count across iterations
# SUPER slow - consider alternative apply version or c++ or parallel
foo <- array(NA, dim = c(nplots, 8, 2))
for(i in 1:nplots) {
  for(j in 1:8) {
    for(k in 1:2) {
      foo[i,j,k] <- mean(unlist(dm100p.MYGA[ , c(paste0("count.hat[", i,",",j,",", k, "]"))]))
    }
  }
}
