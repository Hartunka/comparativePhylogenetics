## run main script to fit pairwise trait models for given family

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

lang.fam.names <- c( "Afro-Asiatic", "Algic", "Arawakan", "Athapaskan-Eyak-Tlingit", "Atlantic-Congo",
                     "Austroasiatic", "Austronesian", "Cariban", "Central_Sudanic", "Chibchan", 
                     "Dravidian", "Gunwinyguan", "Indo-European", "Mande", "Mongolic",  
                     "Nakh-Daghestanian", "Nilotic", "Nuclear_Torricelli", "Nuclear_Trans_New_Guinea", 
                     "Otomanguean", "Pama-Nyungan", "Pano-Tacanan", "Salishan", "Sepik", "Sino-Tibetan",
                     "Siouan", "Sko", "Surmic", "Tai-Kadai", "Tucanoan", "Tupian", 
                     "Turkic", "Uralic",
                     "Uto-Aztecan" ,
                     "Bantu"
                     )

trait.names <- c( "AN", "PN", "ND", "NG", "NNum", "VO", "NRc", "VS" )
t.pairs <- c()

# ---- main loop ----

for (trait.name.1 in trait.names){
  for (trait.name.2 in trait.names){
    pair = paste( trait.name.1, trait.name.2 )
    
    if ( (trait.name.1 != trait.name.2) & !(pair %in% t.pairs) ){
      print( paste("----", trait.name.1, "-", trait.name.2, "----") )
      
      t.pairs = c( t.pairs, pair, paste(trait.name.2, trait.name.1) )
      
      for ( i in 1:length(lang.fam.names) ){
        
        # ---- prepare data ----
          
        lang <- lang.fam.names[i]
        
        print( paste("---", lang, "---") )
        
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
        #dat$C <- avg.Bantu
        
        # ---- prepare metrics ----
        
        save.dir.looc <- file.path(post.dir, "looc", lang)
        save.dir.waic <- file.path(post.dir, "waic", lang)
        save.dir.bf <- file.path(post.dir, "bridge", lang)
        if (!dir.exists(save.dir.looc)) { dir.create(save.dir.looc) }
        if (!dir.exists(save.dir.waic)) { dir.create(save.dir.waic) }
        if (!dir.exists(save.dir.bf)) { dir.create(save.dir.bf) }
        
        # ============== dependent binomial model ==============
          
        print( "-- dependent --" )
        file.name <- file.path( paste(model.name, "-dep_", data.name, "-0", sep="") )
        
        fit.rstan( save.dir=save.dir, file.name=file.name, model=model.dep, dat=dat, 
                   loadExisting=FALSE, seed=1,
                   chains=chains, iter=iter, refresh=2000,
                   multi.C="_5" )
        print("...fit")
        gc(verbose=FALSE)
        
        # ---- eval ----
        get.bridge( model=model.dep, lang.name=lang, 
                    trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                    variant="-dep_", dat=dat, cores=chains,
                    multi.C="_5" )
        print("...bridge")
        gc(verbose=FALSE)
        get.LOOCs( lang.name=lang, 
                   trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                   variant="-dep_", cores=chains,
                   multi.C="_5" )
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
        gc(verbose=FALSE)
        
        # ---- eval ----
        get.bridge( model=model.indep, lang.name=lang, 
                    trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                    variant="-indep_", dat=dat, cores=chains,
                    multi.C="_5" )
        print("...bridge")
        gc(verbose=FALSE)
        get.LOOCs( lang.name=lang, 
                   trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                   variant="-indep_", cores=chains,
                   multi.C="_5" )
        print("...looc")
        gc(verbose=FALSE)
        
        remove( dat, lang.C )
        gc(verbose=FALSE)
      }
    }
  }
}