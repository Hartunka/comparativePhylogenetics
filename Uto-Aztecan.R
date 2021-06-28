## run main script to fit pairwise trait models for Uto-Aztecan languages

# load functions ----
source('paths.R', encoding="utf8" )
source( file.path(fit.dir, 'fit.R'), encoding="utf8" )
source( file.path(dataprep.dir, 'dataprep.R'), encoding="utf8" )
source( file.path(post.dir, 'diagnostics.R'), encoding="utf8" )

library(rstan)
library(bridgesampling)

# prepare models and metaparameters ----
chains <- 14
#max.C <- 100
iter <- 1500

print("parse models...")
model.dep <- rstan::stan_model( 
  file.path( models.dir, 'categorical', 'TwoTraits-binom-multiC-cholesky-nonCentered-main.stan' ),
  #file.path( models.dir, 'categorical', 'Binom-cholesky-nonCentered-main.stan' ),
  model_name="cholesky-binom_dep" )
print("...dep done")
model.indep <- rstan::stan_model( 
  file.path( models.dir, 'categorical', 'TwoTraits-indep-binom-multiC-cholesky-nonCentered-main.stan' ),
  #file.path( models.dir, 'categorical', 'Binom-cholesky-nonCentered-indep-main.stan' ),
  model_name="cholesky-binom_indep" )
print("...indep done")

lang.fam.names <- c( "Uto-Aztecan" )
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
      
      for ( i in 1:length(lang.fam.names) ){
        
        # prepare data for current pair ----
          
        lang <- lang.fam.names[i]
        
        # print( paste("--", lang) )
        
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
        dat$C <- lang.C#[[1]]
        dat$M <- length(lang.C)
          
        # ============== dependent binomial model ==============
          
        print( "- dependent" )
        file.name <- file.path( paste(model.name, "-dep_", data.name, "-0", sep="") )
        
        fit.rstan( save.dir=save.dir, file.name=file.name, model=model.dep, dat=dat, 
                   loadExisting=FALSE, seed=1,
                   chains=chains, iter=iter, refresh=2000,
                   multi.C="_5" )
        print("...fit")
        # ---- eval ----
        get.bridge( model=model.dep, lang.name=lang, 
                    trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                    variant="-dep_", dat=dat, cores=chains, multi.C="_5" )
        print("...bridge")
        gc(verbose=FALSE)
        get.LOOCs( lang.name=lang, 
                   trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                   variant="-dep_", cores=chains, multi.C="_5" )
        print("...looc")
        gc(verbose=FALSE)
        
        # ============== independent binomial model ==============
        
        print("")
        print( "- independent" )
        file.name <- file.path( paste(model.name, "-indep_", data.name, "-0", sep="") )
        
        fit.rstan( save.dir=save.dir, file.name=file.name, model=model.indep, dat=dat, 
                   loadExisting=FALSE, seed=1,
                   chains=chains, iter=iter, refresh=2000,
                   multi.C="_5" )
        print("...fit")
        # ---- eval ----
        get.bridge( model=model.indep, lang.name=lang, 
                    trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                    variant="-indep_", dat=dat, cores=chains, multi.C="_5" )
        print("...bridge")
        gc(verbose=FALSE)
        get.LOOCs( lang.name=lang, 
                   trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                   variant="-indep_", cores=chains, multi.C="_5" )
        print("...looc")
        gc(verbose=FALSE)
          
        remove( dat, lang.C )
        gc(verbose=FALSE)
      }
    }
  }
}