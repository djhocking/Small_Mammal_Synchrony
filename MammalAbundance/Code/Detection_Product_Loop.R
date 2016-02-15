http://stackoverflow.com/questions/20134031/how-to-write-for-loop-when-function-increases-with-each-iteration

set.seed(123)
testp <- matrix(runif(108, 0.1, 0.5), 108, 5)
testmu <- matrix(NA, 108, 5)

for(i in 1:nsites){
  testmu[i,1] <- testp[i,1]
  testmu[i,2] <- testp[i,2]*(1-testp[i,1])
  testmu[i,3] <- testp[i,3]*(1-testp[i,1])*(1-testp[i,2])
  testmu[i,4] <- testp[i,4]*(1-testp[i,1])*(1-testp[i,2])*(1-testp[i,3])
  testmu[i,5] <- testp[i,5]*(1-testp[i,1])*(1-testp[i,2])*(1-testp[i,3])*(1-testp[i,4])
}

testmu2 <- testp*t(apply(cbind(1,1-testp[,-5]),1,cumprod))


testmu2 <- matrix(NA, 108, 5)
nsites = 108
np = 5

# BUGS solution
for (i in 1:nsites) {
  fac <- 1
  testmu2[i,1] <- testp[i,1]
  for (j in 2:np) {
    fac <- fac * (1-testp[i,j-1])
    testmu2[i,j] <- testp[i,j] * fac
  }
}



for (i in 1:nsites) {
  fac <- 1
  testmu2[i,1] <- testp[i,1]
  for (j in 2:np) {
    testmu2[i,j] <- testp[i,j]*(fac <- fac * (1-testp[i,j-1]))
  }
}

testmu3 <- matrix(NA, 108, 5)
nsites = 108
np = 5

for (i in 1:nsites) {
  testmu3[ i, ] <- Reduce( function(x,y) x*(1-y), testp[i, ], 
                           accumulate=TRUE)
}

testmu3



t(apply(testp,1,function(z) Reduce( function(x,y) x*(1-y),z,accumulate=TRUE)))