## run waic comparison(s) for specified stanfit models

# load functions ----
source('paths.R', encoding="utf8" )
source( file.path(fit.dir, 'fit.R'), encoding="utf8" )
source( file.path(dataprep.dir, 'dataprep.R'), encoding="utf8" )
source( file.path(post.dir, 'diagnostics.R'), encoding="utf8" )

library(rstan)
library(rethinking)

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
      
      model.name <- "cholesky-binom"
      data.name <- paste( lang, "_" ,trait.name.1, "-", trait.name.2, sep="" )
        
      file.name <- file.path( paste(model.name, "-lin_", data.name, "-0", sep="") )
      
      waic <- get.WAICs( model.1.name="cholesky-binom-lin",
                         model.2.name="cholesky-binom-univ",
                         lang.name=lang, 
                         trait.name.1=trait.name.1, trait.name.2=trait.name.2,
                         multi.C="_5")
      gc(verbose=FALSE)
      print("...lin/univ")
      
      waic <- get.WAICs( model.1.name="cholesky-binom-univ",
                         model.2.name="cholesky-binom-univ-indep",
                         lang.name=lang, dep="dep_",
                         trait.name.1=trait.name.1, trait.name.2=trait.name.2,
                         multi.C="_5")
      gc(verbose=FALSE)
      print("...univ dep/indep")
      
      gc(verbose=FALSE)
      print("...three")
      
    }
  }
}