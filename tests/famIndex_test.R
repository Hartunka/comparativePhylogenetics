# read 'sub'matrices from covariance matrix over all languages
dat.N <- c(58, 5, 13, 8, 114, 20, 152, 5, 19, 6, 7, 5, 59, 11, 7, 12, 
           14, 7, 31, 17, 22, 6, 5, 6, 74, 7, 5, 8, 6, 7, 11, 10, 8, 23)


ranges <- matrix(dat=0, nrow=34, ncol=4)
ranges[1,1] <- 1;
ranges[1,3] <- 1;
ranges[1,2] <- dat.N[1];
ranges[1,4] <- (2*dat.N[1]);
for (f in 2:34){
  ranges[f,1] <- ranges[f-1,2]+1;
  ranges[f,2] <- ranges[f,1]+dat.N[f]-1;
  ranges[f,3] <- ranges[f-1,4]+1;
  ranges[f,4] <- ranges[f,3]+(2*dat.N[f])-1;
}

p <- vector(mode="character", length=(2*768))
Cmat <- matrix(data=0, nrow=768, ncol=768)
for (i in 1:768) {
  for (j in 1:768) {
    Cmat[i,j] <- paste(i,j,sep=",")
  }
}
Cs <- vector(mode="list", length=34)

for (f in 1:34){
  p[ranges[f,3]:ranges[f,4]] <- paste(f,": ",dat.N[f],"(",dat.N[f]*2,"); ", 
                                     ranges[f,3], ":", ranges[f,4],
                                     sep="");
  C.f <- Cmat[ ranges[f,1]:ranges[f,2], ranges[f,1]:ranges[f,2] ]
  Cs[[f]] <- C.f
}