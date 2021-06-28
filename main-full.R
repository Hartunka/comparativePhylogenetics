## run main script to fit pairwise trait models over all languages

# load functions ----
source('paths.R', encoding="utf8" )
source( file.path(fit.dir, 'fit.R'), encoding="utf8" )
source( file.path(dataprep.dir, 'dataprep.R'), encoding="utf8" )
source( file.path(post.dir, 'diagnostics.R'), encoding="utf8" )

library(rstan)
library(bridgesampling)

# prepare models and metaparameters ----
chains <- 14
iter <- 1500

print("parse models...")
model.lin <- rstan::stan_model( 
  file.path( models.dir, 'categorical', 'Binom-cholesky-nonCentered-fullLang-multiC-lin.stan' ),
  #file.path( models.dir, 'categorical', 'Binom-cholesky-nonCentered-fullLang-lin.stan' ),
  model_name="cholesky-binom_lin" )
print("...lineage done")
model.univ <- rstan::stan_model( 
  file.path( models.dir, 'categorical', 'Binom-cholesky-nonCentered-fullLang-multiC-univ.stan' ),
  #file.path( models.dir, 'categorical', 'Binom-cholesky-nonCentered-fullLang-univ.stan' ),
  model_name="cholesky-binom_univ" )
print("...universal dep done")
#model.univ.indep <- rstan::stan_model( 
#  file.path( models.dir, 'categorical', 'Binom-cholesky-nonCentered-fullLang-univ-indep.stan' ),
#  model_name="cholesky-binom_univ-indep" )
#print("...universal indep done")

trait.names.a <- c( "AN", "PN", "ND", "NG", "NNum", "VO", "NRc", "VS" )
trait.names.b <- c( "PN", "ND", "NG", "NNum", "VO", "NRc", "VS" )
t.pairs <- c()

# ---- main loop ----

for (trait.name.1 in trait.names.a){
  for (trait.name.2 in trait.names.b){
    pair = paste( trait.name.1, trait.name.2 )
    
    if ( (trait.name.1 != trait.name.2) & !(pair %in% t.pairs) ){
      print( paste("---", trait.name.1, "-", trait.name.2) )
      
      t.pairs = c( t.pairs, pair, paste(trait.name.2, trait.name.1) )
      
      # prepare data for current pair ----
          
      lang <- "full"
      
      # setup save
      save.dir <- file.path("fit", "fit", lang)
      if (!dir.exists(save.dir)) { dir.create(save.dir) }
      model.name <- "cholesky-binom"
      data.name <- paste( lang, "_" ,trait.name.1, "-", trait.name.2, sep="" )
        
      # load prepared data as "dat"
      data.file.name <- file.path( dataprep.dir, "dat", lang, 
                                   paste(data.name, "binom-0_dat.RData", sep="_") )
      C.file.name <- file.path( dataprep.dir, "dat", lang, 
                                   paste(lang, "_avgC-5.RData", sep="_") )
      load(data.file.name)
      load(C.file.name)
      dat$C <- lang.C
      
      # ============== lineage specific binomial model ==============
        
      print( "- lineage" )
      file.name <- file.path( paste(model.name, "-lin_", data.name, "-0", sep="") )
        
      fit.rstan( save.dir=save.dir, file.name=file.name, model=model.lin, dat=dat, 
                 loadExisting=FALSE, seed=1,
                 chains=chains, iter=iter, refresh=2000,
                 multi.C="_5" )
      print("...fit")
      get.bridge( model=model.lin, lang.name=lang, 
                  trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                  variant="-lin_", dat=dat, cores=chains,
                  multi.C="_5" )
      gc(verbose=FALSE)
      print("...bridge")
      get.LOOCs( lang.name=lang, 
                 trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                 variant="-lin_", cores=chains,
                 multi.C="_5" )
      gc(verbose=FALSE)
      print("...looc")
          
      # ============== universal-dependent binomial model ==============
          
      print("")
      print( "- universal - dep" )
      file.name <- file.path( paste(model.name, "-univ_", data.name, "-0", sep="") )
        
      fit.rstan( save.dir=save.dir, file.name=file.name, model=model.univ, dat=dat, 
                 loadExisting=FALSE, seed=1,
                 chains=chains, iter=iter, refresh=2000,
                 multi.C="_5" )
      print("...fit")
      get.bridge( model=model.univ, lang.name=lang, 
                  trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                  variant="-univ_", dat=dat, cores=chains,
                  multi.C="_5" )
      print("...bridge")
      gc(verbose=FALSE)
      get.LOOCs( lang.name=lang, 
                 trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                 variant="-univ_", cores=chains,
                 multi.C="_5" )
      gc(verbose=FALSE)
      print("...looc")
      
      # ============== universal-independent binomial model ==============
      
      #print("")
      #print( "- universal - indep" )
      #file.name <- file.path( paste(model.name, "-univ-indep_", data.name, "-0", sep="") )
      
      #fit.rstan( save.dir=save.dir, file.name=file.name, model=model.univ.indep, dat=dat, 
      #           loadExisting=FALSE, seed=1,
      #           chains=chains, iter=iter, refresh=2000 )
      #print("...fit")
      
      remove( dat, lang.C )
      gc(verbose=FALSE)
    }
  }
}