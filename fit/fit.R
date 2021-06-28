
#' create design matrix to be used for Multivariate normal distribution's mean
#' 
#' @param N number of languages
#' @param r number of feats
#' @return dMat design Matrix
#' @export
#'
designMatrix <- function(N, r){
  dMat <- matrix(0, nrow=N*r, ncol=r)
  for (i in 1:(N*r)) {
    for (j in 1:r) {
      if ( (j-1*N < i)&(i<=j*N) ){
        dMat[i,j] <- 1
      }
    }
  }
  return(dMat)
}

## ---- main ----

#' fit specified data to stan model via rstan
#' 
#' if fit model object already exists at specified save location, 
#' don't run stan and load existing object instead
#' 
#' 
#' @param save.dir character, path to write stanfit object to
#' @param file.name file name for stanfit object
#' @param model compiled stanmodel object to be fit
#' @param dat data to be fit to, as list
#' @param loadExisting if TRUE and stanfit already exists, load into workspace
#' @param seed random seed for stan's sampling
#' @param chains number of MCMC chains to run
#' @param refresh iteration interval at for which to print status update
#' @param iter number of iterations after warmup
#' @param adapt_delta,max_treedepth  parameters to be passed to stan
#' @param multiC covariance matrix variant used as string for naming resulting file
#' @return list:
#'         fit: stanfit object
#'         t.diff: difftime object, result of taking time for model run 
#'                 (or file loading)
#'         path: location, where stanfit object has been written to
#'         file.name: same as input file.name (for convenience)
#'
fit.rstan <- function( save.dir, file.name, model, dat, loadExisting=TRUE,
                       seed=1, chains=4, refresh=1600, iter=4000,
                       adapt_delta=0.8, max_treedepth=10, multi.C="" ){
  
  save.path.fit <- file.path(save.dir, paste(file.name, multi.C, "_stanfit.RData", sep=""))
  print(save.path.fit)
  start.time <- Sys.time()
  if ( !file.exists(save.path.fit) ){
    
    stanfit <- rstan::sampling(object=model, data=dat, 
                       iter=(iter+1500), warmup=1500, chain=chains, cores=chains,
                       refresh=refresh, seed=seed,
                       save_warmup=TRUE,
                       pars=c("V"), include=FALSE)
    # save model fit
    save( stanfit, file=save.path.fit )
  } else {
    if (loadExisting){
      print("already exists, loading stanfit...")
      load( file=save.path.fit )
      
    } else{
      print("already exists, return")
      stanfit<-NA
    }
  }
  end.time <- Sys.time()
  print(end.time - start.time)
  t.diff <- difftime( end.time, start.time, units=c("mins") )
  
  return( list(fit=stanfit, 
               t.diff=t.diff, path=save.path.fit, name=file.name ) )
}
## ---- cross validation ----
# eventually not used

fit.CV.2f <- function( model, dat.1, dat.2, indep=FALSE, 
                       lang, trait.name.1, trait.name.2 ){
  
  if (indep){
    depStr <- "_indep"
  } else {
    depStr <- ""
  }
  
  # -- fold 1 -- 
  mod.1.name <- paste(lang, "_", trait.name.1, "-", trait.name.2, "_", dat.1$M, depStr, "_CV-f1", sep="")
  stanfit.1.path <- file.path( "CV-out", lang, paste(mod.1.name, "RData", sep=".") )
  
  if (file.exists(stanfit.1.path)){
    load( file=stanfit.1.path )
  } else {
    
    fit.1 <- model$sample(
      data=dat.1, chains=6, parallel_chains=6, iter_sampling=2000, iter_warmup=2000,
      refresh=2000, save_warmup=FALSE
    )
    
    suppressMessages( expr=fit.1$save_output_files( dir="CV-out", basename=mod.1.name, 
                                                    timestamp=FALSE, random=FALSE ) 
    )
    stanfit.1 <- rstan::read_stan_csv( fit.1$output_files() )
    save( stanfit.1, file=stanfit.1.path )
    
    rm(fit.1)
    suppressMessages( gc(verbose=FALSE) )
  }
  
  # get accuracy
  samples.1 <- rstan::extract(stanfit.1, pars=c("accuracy"))
  
  acc.1 <- mean( samples.1$accuracy ) 
  
  rm(stanfit.1)
  suppressMessages( gc(verbose=FALSE) )
  
  # -- fold 2 -- 
  mod.2.name <- paste(lang, "_", trait.name.1, "-", trait.name.2, "_", dat.2$M, depStr, "_CV-f2", sep="")
  stanfit.2.path <- file.path( "CV-out", lang, paste(mod.2.name, "RData", sep=".") )
  
  if (file.exists(stanfit.2.path)){
    load( file=stanfit.2.path )
  } else {
    
    fit.2 <- model$sample(
      data=dat.2, chains=6, parallel_chains=6, iter_sampling=2000, iter_warmup=2000,
      refresh=2000, save_warmup=FALSE
    )
    
    suppressMessages( expr=fit.2$save_output_files( dir="CV-out", basename=mod.2.name, 
                                                    timestamp=FALSE, random=FALSE ) 
    )
    stanfit.2 <- rstan::read_stan_csv( fit.2$output_files() )
    save( stanfit.2, file=stanfit.2.path )
    
    rm(fit.2)
    suppressMessages( gc(verbose=FALSE) )
  }
  
  # get accuracy
  samples.2 <- rstan::extract(stanfit.2, pars=c("accuracy"))
  
  acc.2 <- mean( samples.2$accuracy ) 
  
  rm(stanfit.2)
  suppressMessages( gc(verbose=FALSE) )
  
  return( list(acc.1=acc.1, acc.2=acc.2, acc=mean(acc.1, acc.2)) )
}


fit.CV.3f <- function( model, dat.1, dat.2, dat.3, indep=FALSE, 
                       lang, trait.name.1, trait.name.2 ){
  
  if (indep){
    depStr <- "_indep"
  } else {
    depStr <- ""
  }
  
  # -- fold 1 -- 
  mod.1.name <- paste(lang, "_", trait.name.1, "-", trait.name.2, "_", dat.1$M, depStr, "_CV-f1", sep="")
  stanfit.1.path <- file.path( "CV-out", lang, paste(mod.1.name, "RData", sep=".") )
  
  if (file.exists(stanfit.1.path)){
    load( file=stanfit.1.path )
  } else {
    
    fit.1 <- model$sample(
      data=dat.1, chains=6, parallel_chains=6, iter_sampling=2000, iter_warmup=2000,
      refresh=2000, save_warmup=FALSE
    )
    
    suppressMessages( expr=fit.1$save_output_files( dir="CV-out", basename=mod.1.name, 
                                                    timestamp=FALSE, random=FALSE ) 
    )
    stanfit.1 <- rstan::read_stan_csv( fit.1$output_files() )
    save( stanfit.1, file=stanfit.1.path )
    
    rm(fit.1)
    suppressMessages( gc(verbose=FALSE) )
  }
  
  # get accuracy
  samples.1 <- rstan::extract(stanfit.1, pars=c("accuracy"))
  
  acc.1 <- mean( samples.1$accuracy ) 
  
  rm(stanfit.1)
  suppressMessages( gc(verbose=FALSE) )
  
  # -- fold 2 -- 
  mod.2.name <- paste(lang, "_", trait.name.1, "-", trait.name.2, "_", dat.2$M, depStr, "_CV-f2", sep="")
  stanfit.2.path <- file.path( "CV-out", lang, paste(mod.2.name, "RData", sep=".") )
  
  if (file.exists(stanfit.2.path)){
    load( file=stanfit.2.path )
  } else {
    
    fit.2 <- model$sample(
      data=dat.2, chains=6, parallel_chains=6, iter_sampling=2000, iter_warmup=2000,
      refresh=2000, save_warmup=FALSE
    )
    
    suppressMessages( expr=fit.2$save_output_files( dir="CV-out", basename=mod.2.name, 
                                                    timestamp=FALSE, random=FALSE ) 
    )
    stanfit.2 <- rstan::read_stan_csv( fit.2$output_files() )
    save( stanfit.2, file=stanfit.2.path )
    
    rm(fit.2)
    suppressMessages( gc(verbose=FALSE) )
  }
  
  # get accuracy
  samples.2 <- rstan::extract(stanfit.2, pars=c("accuracy"))
  
  acc.2 <- mean( samples.2$accuracy ) 
  
  rm(stanfit.2)
  suppressMessages( gc(verbose=FALSE) )
  
  # -- fold 3 -- 
  mod.3.name <- paste(lang, "_", trait.name.1, "-", trait.name.2, "_", dat.3$M, depStr, "_CV-f3", sep="")
  stanfit.3.path <- file.path( "CV-out", lang, paste(mod.3.name, "RData", sep=".") )
  
  if (file.exists(stanfit.3.path)){
    load( file=stanfit.3.path )
  } else {
    
    fit.3 <- model$sample(
      data=dat.3, chains=6, parallel_chains=6, iter_sampling=2000, iter_warmup=2000,
      refresh=2000, save_warmup=FALSE
    )
    
    suppressMessages( expr=fit.3$save_output_files( dir="CV-out", basename=mod.3.name, 
                                                    timestamp=FALSE, random=FALSE ) 
    )
    stanfit.3 <- rstan::read_stan_csv( fit.3$output_files() )
    save( stanfit.3, file=stanfit.3.path )
    
    rm(fit.3)
    suppressMessages( gc(verbose=FALSE) )
  }
  
  # get accuracy
  samples.3 <- rstan::extract(stanfit.3, pars=c("accuracy"))
  
  acc.3 <- mean( samples.3$accuracy ) 
  
  rm(stanfit.3)
  suppressMessages( gc(verbose=FALSE) )
  
  return( list(acc.1=acc.1, acc.2=acc.2, acc.3=acc.3, acc=mean(acc.1, acc.2, acc.3)) )
}
