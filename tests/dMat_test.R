# implement design matrix to be used for multinormal mean
# 'D is an nr × r design matrix where each element D_ij is 1 if (j−1)·n < i · j · n and 0 otherwise'
#   Harmon, 2019

n <- 6
r <- 2

dMat <- matrix(0, nrow=n*r, ncol=r)
for (i in 1:(n*r)) {
  for (j in 1:r) {
    if ( (j-1*n < i)&(i<=j*n) ){
      dMat[i,j] <- 1
    }
  }
}
dMat