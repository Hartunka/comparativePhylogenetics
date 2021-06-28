

#' check if specified directory exists and if not, create it
#' @param path to directory to be checked (character)
#' 
check.directory <- function( path ){
  if ( !dir.exists( path ) ){
    dir.create( path )
  }
}

# load stanfit object from specifed location
load.stanfit <- function(save.dir, file.name, cVar=""){
  load( file=file.path(save.dir, paste(file.name, cVar, "_stanfit.RData", sep="")) )
  return(stanfit)
}

#' perform bridgesampling for specified model and save result
#' if file for results already exists, skip sampling and load file instead
#' 
#' @param model stanmodel object to be used for prior evaluation
#' @param lang.name name of language family for which current model has been fit
#' @param trait.name.1 name of first trait for which model has been fit
#' @param trait.name 2 name of second trait for which model has been fit
#' @param variant string identifier for which model variant has been used (dep/indep/univ etc.)
#' @param dat data list for empty model fit for prior evaluation
#' @param cores number of cores to be used for bridgesampling
#' @param multi.C string identifier for which phylogenetic model variant has been used
#' @param trait.var string identifier for whether 'other' category in trait data has been subsumed to 0 or 1
#' 
#' @return bridge object containing results from bridgesampling
#' 
get.bridge <- function( model, lang.name, trait.name.1, trait.name.2, variant, dat, cores=1, multi.C="", trait.var="-0" ){
  
  # load fit model
  save.dir <- file.path( fit.dir, "fit", lang )
  data.name <- paste( lang.name, "_" ,trait.name.1, "-", trait.name.2, sep="" )
  file.name <- file.path( paste("cholesky-binom", variant, data.name, trait.var, sep="") )
  
  if(file.exists(file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")))){
    print( file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")) )
    load( file=file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")) )
    
  } else {
    print( file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")) )
    stanfit <- load.stanfit( save.dir, file.name, cVar=multi.C  )
    
    # create empty stanfit.model
    cat("load empty model...\n")
    stan.model.fit <- sampling( object=model, data=dat, iter=0, chains=1, verbose=FALSE )
    
    # bridgesampling
    cat("bridgesampling...\n")
    bridge <- bridge_sampler( samples=stanfit, stanfit_model=stan.model.fit,
                              silent=TRUE, cores=cores, method="warp3" )
    cat("save...\n")
    save(bridge, file=file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")) )
  }
  return(bridge)
}

#' compute leave-one-out cross validation approximation for a model
#' 
#' @param fit fit model to be evaluated (stanfit)
#' @return fit.loo : list containing evaluations for model returned from \code{loo::loo}
#' 
eval.loo <- function( fit, cores=1 ){
  
  cat( "...loo...\n" )
  
  log_lik <- loo::extract_log_lik(fit, merge_chains=FALSE)
  r_eff <- loo::relative_eff( exp(log_lik), cores=cores )
  fit.loo <- loo::loo( log_lik, r_eff=r_eff, cores=cores )
  
  return( fit.loo )
}

#' get and save LOO object for specified stanfit model
#' if result file already exists, skip computation and load file
#' wrapper for eval.loo()
#' 
#' @param lang.name name of family for which model has been fit
#' @param trait.name.1 name of first trait for which model has been fit
#' @param trait.name 2 name of second trait for which model has been fit
#' @param variant string identifier for which model variant has been used (dep/indep/univ etc.)
#' @param cores number of cores to be used for computation
#' @param multi.C string identifier for which phylogenetic model variant has been used
#' 
#' @return list containing results from loo::loo() applied to specified model
#' 
get.LOOCs <- function( lang.name, trait.name.1, trait.name.2, variant, cores=1, multi.C="" ){
  
  # load fit model
  save.dir <- file.path( fit.dir, "fit", lang )
  data.name <- paste( lang.name, "_" ,trait.name.1, "-", trait.name.2, sep="" )
  file.name <- file.path( paste( "cholesky-binom", variant, data.name, "-0", multi.C, sep="" ) )
  
  cat(file.name, "\n")
  if(file.exists(file.path(post.dir, "looc", lang.name, paste(file.name, "_looc.RData", sep="")))){
    
    load( file=file.path(post.dir, "looc", lang.name, paste(file.name, "_looc.RData", sep="")) )
    
  } else {
    
    stanfit <- load.stanfit(save.dir, file.name, cVar=multi.C )
    # LOOC
    loo <- eval.loo( stanfit, cores=cores )
    save( loo, file=file.path( post.dir, "looc", lang.name, paste(file.name, multi.C, "_looc.RData", sep="") ) )
  }
  return(loo)
}


#' get and save WAIC object from comparing specified stanfit models
#' if result file already exists, skip computation and load file
#' 
#' @param model.1.name of first stanfit model
#' @param model.2.name of second stanfit model
#' @param lang.name name of family for which models have been fit
#' @param trait.name.1 name of first trait for which models have been fit
#' @param trait.name 2 name of second trait for which models have been fit
#' @param multi.C string identifier for which phylogenetic model variant has been used
#' @param dep additional string identifier to discriminate universal/lineage and dependent/independent
#'            comparisons over all languages when saving results
#' 
#' @return compareIC object resulting from applying rethinking::compare() to specified models
#' 
get.WAICs <- function( model.1.name, model.2.name, lang.name, trait.name.1, trait.name.2, multi.C="", dep="_" ){
  
  # load fit model
  save.dir <- file.path( fit.dir, "fit", lang )
  data.name <- paste( lang.name, "_" ,trait.name.1, "-", trait.name.2, sep="" )
  file.name.waic <- file.path( paste( "cholesky-binom", dep, data.name, "-0", multi.C, "_waic.RData", sep="" ) )
  
  cat(file.name.waic, "\n")
  if(file.exists(file.path(post.dir, "waic", lang.name, file.name.waic))){
    
    load( file=file.path(post.dir, "waic", lang.name, file.name.waic) )
    
  } else {
    file.name.1 <- file.path( paste( model.1.name, "_", data.name, "-0", sep="" ) )
    file.name.2 <- file.path( paste( model.2.name, "_", data.name, "-0", sep="" ) )
    
    cat("loading", file.name.1, "...\n")
    stanfit.1 <- load.stanfit(save.dir, file.name.1, cVar=multi.C)
    cat("loading", file.name.2, "...\n")
    stanfit.2 <- load.stanfit(save.dir, file.name.2, cVar=multi.C)
    
    # WAIC
    print("waic...")
    waic <- rethinking::compare( stanfit.1, stanfit.2 )
    
    save( waic, file=file.path( post.dir, "waic", lang.name, file.name.waic ) )
  }
  return( waic )
}


#' plot correlated pairs in the style of Dunn et al., Jäger or in custom style
#' 
#' @param correlations matrix containing correlations between trait pairs
#' @param threshold minimum value for correlations to appear in plot
#' @param as.log plot values on log modified log scale for visual clarity
#' @param log.scale modification for log scale (multiplier)
#' @param comp.corr additional correlation matrix to compare main correlations to
#' @param unif.traits traits to bey greyed in style of Dunn et al.
#' @param fit.colour colour for cases where compared correlations agree
#' @param miss.colour colour for cases where compared correlations disagree
#' @param max.width maximum width of bars presenting correlations
#' @param outline how to arrange traits on plot
#' 
#' @return plot representing correlated pairs
#' 

dunn.plot <- function(correlations, threshold=5, as.log=FALSE, log.scale=5, comp.corr=NA, unif.traits=NA,
                      fit.colour="green", miss.colour="red", max.width=20, outline="dunn"){
  
  `%!in%` <- Negate(`%in%`)
  # -- setup focal points in df --
  if (outline=="dunn"){
    pos.mat <- matrix(data=c(0,2, 0,4, 2,0, 2,6, 
                             5.5,0, 5.5,6, 7.5,2, 7.5,4), 
                      nrow=8, byrow=TRUE )
    trait.df <- as.data.frame(pos.mat)
    colnames(trait.df) <- c("x", "y")
    trait.df$labels <- c("Adjective\n-noun", "Subject\n-verb", "Demonstrative\n-noun", "Numeral\n-noun", 
                         "Relative\nclause-noun", "Genetive\n-noun", "Object\n-verb", "Adposition\n-noun")
    trait.df$l.short <- c("AN", "VS", "ND", "NNum", "NRc", "NG", "VO", "PN")
    
  }
  if (outline=="jaeger"){
    pos.mat <- matrix(data=c(0,0, 6,0, 11,0.5, 
                             4,3, 10,3,
                             3,6, 11,6, 
                             7,8), 
                      nrow=8, byrow=TRUE )
    trait.df <- as.data.frame(pos.mat)
    colnames(trait.df) <- c("x", "y")
    trait.df$labels <- c("Adjective\n-noun", "Demonstrative\n-noun",  "Numeral\n-noun", 
                         "Relative\nclause-noun", "Object\n-verb", 
                         "Subject\n-verb", "Adposition\n-noun",  
                         "Genetive\n-noun")
    trait.df$l.short <- c("AN", "ND", "NNum", "NRc", "VO", "VS", "PN", "NG")
  }
  if (outline=="new"){
    pos.mat <- matrix(data=c(0,0, 0,3, 6,0, 
                             6,3, 9,6,
                             0,9, 3,6, 
                             6,9), 
                      nrow=8, byrow=TRUE )
    trait.df <- as.data.frame(pos.mat)
    colnames(trait.df) <- c("x", "y")
    trait.df$labels <- c("Adjective\n-noun", "Demonstrative\n-noun",  "Numeral\n-noun", 
                         "Relative\nclause-noun", "Object\n-verb", 
                         "Subject\n-verb", "Adposition\n-noun",  
                         "Genetive\n-noun")
    trait.df$l.short <- c("AN", "ND", "NNum", "NRc", "VO", "VS", "PN", "NG")
  }
  
  
  # -- setup basic plot -- 
  corr.plot <- ggplot(trait.df, aes(x=x, y=y))
  
  if (is.matrix(comp.corr)){
    # -- add comparisons --
    is <- c()
    for (i in 1:nrow(trait.df)) {
      point.1 <- c(trait.df[i,1], trait.df[i,2])
      is <- c(is, i)
      for (j in 1:nrow(trait.df)) {
        if (i!=j && !(j %in% is)){
          point.2 <- c(trait.df[j,1], trait.df[j,2])
          
          weight <- correlations[i,j]
          c.weight <- comp.corr[i,j]
          
          if (c.weight>threshold){
            
            c.width <- min(ifelse(as.log, log(c.weight)*log.scale, c.weight),max.width)
            width <- min(ifelse(as.log, log(weight)*log.scale, weight), max.width)
            
            if (weight<threshold){
              corr.plot <- corr.plot + annotate( "segment", 
                                                 x = point.1[1],
                                                 xend = point.2[1],
                                                 y = point.1[2],
                                                 yend = point.2[2],
                                                 size = c.width, 
                                                 alpha = 0.5,
                                                 colour = miss.colour
              )
            } else {
              c.width <- 27
              corr.plot <- corr.plot + annotate( "segment", 
                                                 x = point.1[1],
                                                 xend = point.2[1],
                                                 y = point.1[2],
                                                 yend = point.2[2],
                                                 size = c.width+5,
                                                 alpha = 0.5,
                                                 colour = fit.colour
              )
            }
          }
        }
      }
    }
  }
  
  # -- add correlations --
  is <- c()
  for (i in 1:nrow(trait.df)) {
    point.1 <- c(trait.df[i,1], trait.df[i,2])
    is <- c(is, i)
    for (j in 1:nrow(trait.df)) {
      if (i!=j && !(j %in% is)){
        point.2 <- c(trait.df[j,1], trait.df[j,2])
        
        weight <- correlations[i,j]
        
        if (weight>threshold){
          width <- min(ifelse(as.log, log(weight)*log.scale, weight), max.width)
          corr.plot <- corr.plot + annotate( "segment", 
                                             x = point.1[1],
                                             xend = point.2[1],
                                             y = point.1[2],
                                             yend = point.2[2],
                                             size = width,
                                             alpha = 0.8,
                                             colour = "black"
          )
        }
        
      }
    }
  }
  
  corr.plot <- corr.plot + geom_label( label=trait.df$labels,
                                       label.padding = unit(0.55, "lines"), # Rectangle size around label
                                       label.r = unit(2, "lines"), # Rectangle corner rounding
                                       label.size = 1,
                                       color = "black",
                                       fill = "white",
                                       size = rel(7)
  )
  
  if (is.character(unif.traits)){
    unif.df <- trait.df[trait.df$l.short %in% unif.traits,]
    
    corr.plot <- corr.plot + geom_label( mapping=aes(x=x, y=y),
                                         label=unif.df$labels,
                                         data=unif.df,
                                         label.padding = unit(0.55, "lines"),
                                         label.r = unit(2, "lines"),
                                         label.size = 1,
                                         color = "darkgrey",
                                         fill = "white",
                                         size = rel(7)
    )
  }
  
  corr.plot <- corr.plot + theme( panel.background=element_blank(),
                                  axis.line=element_blank(), axis.ticks=element_blank(),
                                  axis.text.x=element_blank(), axis.text.y=element_blank(),
                                  axis.title.x=element_blank(), axis.title.y=element_blank()
  )
  if (outline=="dunn"){
    corr.plot <- corr.plot + xlim(-0.8, 8.4) + ylim(-0.5, 6.5)
  }
  if (outline=="jaeger"){
    corr.plot <- corr.plot + xlim(-0.8, 12.4) + ylim(-0.5, 8.5)
  } 
  if (outline=="new"){
    corr.plot <- corr.plot + xlim(-2, 10) + ylim(-0.5, 9.5)
  } 
  
  return(corr.plot)
}

