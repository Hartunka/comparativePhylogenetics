
tree.glot <- list(
  "Afro-Asiatic", "Algic", "Arawakan", "Athapaskan-Eyak-Tlingit", "Atlantic-Congo", "Austroasiatic", "Austronesian",
  "Cariban", "Central_Sudanic", "Chibchan",
  "Dravidian", "Gunwinyguan", "Indo-European", 
  "Mande", "Mongolic",
  "Nakh-Daghestanian", "Nilotic", "Nuclear_Torricelli", "Nuclear_Trans_New_Guinea",
  "Otomanguean",
  "Pama-Nyungan", "Pano-Tacanan",
  "Salishan", "Sepik", "Sino-Tibetan", "Siouan", "Sko", "Surmic",
  "Tai-Kadai", "Tucanoan", "Tupian", "Turkic",
  "Uralic", "Uto-Aztecan"
  )

# ======== read features from csv ========

#' read language features from csv file "charMtx.csc"
#' 
#' filters for features of one specified language family 
#' or for all families for which relation information is given in list tree.glot
#' 
#' transforms value types to numeric - if bin==TRUE, binary - values, such that
#' | "-" | 1 | 0/1 |
#' | "0" | 2 | 0   |
#' | "1" | 3 | 1   |
#' 
#' @param lang.name the name of language family to get features for, or "all" for unfiltered
#' @param bin if TRUE features are to be read as binary
#' @param other if bin==TRUE, is the 'other' category "-" to be read as 0 or 1
#' @return feats.lang : a data.frame containing:
#' | X    | family.subfam.lang name    |                          |                          |             |
#' | AN   | noun-adjective order       | "0" : noun-adjective     | "1" : adjective-noun     | "-" : other |
#' | PN   | adposition-noun order      | "0" : prepositions       | "1" : adpositions        | "-" : other |
#' | ND   | noun-demonstrative order   | "0" : noun-demonstrative | "1" : demonstrative-noun | "-" : other |
#' | NG   | noun-genetive order        | "0" : noun-genetive      | "1" : genetive-noun      | "-" : other |
#' | NNum | noun-numeral order         | "0" : noun-numeral       | "1" : numeral-noun       | "-" : other |
#' | VO   | verb-object order          | "0" : verb-object        | "1" : object-verb        | "-" : other |
#' | NRc  | noun-relative clause order | "0" : noun-relative      | "1" : relative-noun      | "-" : other |
#' | VS   | verb-subject order         | "0" : verb-subject       | "1" : subject-verb       | "-" : other |
#' | glot | glottolog lang family name |                          |                          |             |
#' 
read.feats <- function( lang.name, bin=FALSE, other=0 ){
  
  feats <- read.csv( file.path( 
    data.dir,     # path to directory containing data, defined in "paths.R"
    "charMtx.csv" # name of file containing the feature information
    ) )
  # renaming feature name to avoid naming conflict with NA
  names(feats)[names(feats) == 'NA.'] <- 'AN'
  
  if (lang.name=="all"){
    # filter all languages, for which no relation information is given
    feats.filtered <- feats[feats$glot %in% tree.glot,]
    # sort by language family
    order.glot <- order(feats.filtered$glot)
    feats.lang <- feats.filtered[order.glot, ]
  } else {
  # filter specified language
    feats.lang <- feats[feats$glot == lang.name,]
  }
  
  # features as numerical, "-" -> 1, "0" -> 2, "1" -> 3
  # whether they're read as strings or as factor seems to depend on the R version
  if (is.factor(feats.lang$AN)){  # =< R 3.5
    for (column in 2:( ncol(feats.lang)-1) ) {
      feats.lang[column] <- as.numeric( feats.lang[[column]] )
    }
  } else if (is.character(feats.lang$AN)){ # >= R 4.0
    for (column in 2:( ncol(feats.lang)-1) ) {
      # transform first to factor, to fit with how csv is read in older version
      feats.lang[column] <- factor( feats.lang[[column]], c("-", "0", "1") )
      feats.lang[column] <- as.numeric( feats.lang[[column]] )
    }
  }
  # features as binary, "-"/1 -> 0|1, "0"/2 -> 0, "1"/3 -> 1
  if (bin){
    for (column in 2:( ncol(feats.lang)-1) ) {
      feats.lang[[column]][ feats.lang[[column]] != 3-other] <- 0+other
      feats.lang[[column]][ feats.lang[[column]] == 3-other] <- 1-other
    }
  }
  
  return(feats.lang)
}


# ======== read trees ========

#' read phylogenetic tree data for specified language family 
#' from file in data directory specified in paths.R
#'
#' @param lang.name name of language family
#' @param feats.lang.names array of names of languages of family specified with lang.name
#'                         phylogenetic information for languages not in feats.lang.names will be filtered out
#' @param phyl if TRUE, return trees of class phylo, else covariance matrices
#' @param avg if trues computes matrices representing average phylogenetic distances for the specified family
#' @param n_avg number of distinct average matrices to compute if avg==TRUE
#' @return if phyl=FALSE trees.C : list of covariance matrices 
#'                                 for specified language family derived from phylogenetic trees
#'         else trees.trimmed : list of phylo, containing phylogenetic trees for family
#'
read.trees <- function(lang.name, feats.lang.names, 
                       avg=FALSE, n_avg=1, phyl=FALSE ){
                       
  trees <- ape::read.tree(
    file.path( data.dir, "trees", paste(lang.name, ".posterior.tree", sep="" ) )
    )
  # filter trees to only contain languages for which feature information is given
  trees.trimmed <- lapply( trees, 
                           function(i) { ape::keep.tip( i, feats.lang.names ) }
                           )
  if (phyl){
    return(trees.trimmed)
  }
  if (!avg) {
    trees.C <- list()
    for ( i in 1:length( trees.trimmed ) ){
      # suppressWarnings because otherwise corBownian keeps returning 
      # messages about using the default value for 'form'parameter
      suppressWarnings({
        tree <- trees.trimmed[[i]] # sort matrices, so each matrix lists languages in the same order:
        tree.C <- ape::vcv( ape::corBrownian(1, tree), corr=FALSE )[feats.lang.names, feats.lang.names]
        trees.C[[i]] <- tree.C
      })
    }
  return( trees.C )
  } else{
    nC <- length(trees.trimmed)
    n.split <- ceiling(nC/n_avg)
    tree.splits <- split(trees.trimmed, rep(1:ceiling(nC/n.split), each=n.split, length.out=nC))
    avg.Cs <- list()
    
    for ( h in 1:length(tree.splits) ) {
      tree.split <- tree.splits[[h]]
      
      nrc <- length(tree.split[[1]]$tip.label)
      sum.C <- matrix(data=rep(0, nrc*nrc),  nrow=nrc, ncol=nrc)
      for ( i in 1:length( tree.split ) ) {
        suppressWarnings({
          tree <- tree.split[[i]]
          C.i <- ape::vcv( ape::corBrownian(1, tree), corr=FALSE )[feats.lang.names, feats.lang.names]
        })
        sum.C = sum.C + C.i
      }
      avg.Cs[[h]] = sum.C / length( tree.split )
    }
    return(avg.Cs)
  }
}

#' create covariance matrix for all language families, 
#' or a list of languages from input
#' 
#' if avg==TRUE for each language family the values are the average over all trees for the family
#' 
#' @param feats data frame as read from read.feats() to provide information about feature order and family names
#' @param avg if TRUE return average covariance matrix
#' @param n_avg number of distinct average matrices to compute if n_avg==TRUE
#' @return (average) covariance matrix for all or specified languages families
#' 
read.trees.multi.Lang <- function(feats, avg=TRUE, n_avg=1){
  
  fam.names <- unique(feats$glot)
  trees <- list()
  for (i in 1:length(fam.names)) {
    lang.names <- feats[feats$glot==fam.names[i],]$X
    trees_i <- read.trees(fam.names[i], lang.names, avg=avg, n_avg=n_avg)
    if (!avg){
      set.seed(1)
      trees[[i]] <- sample(trees_i, max.C)
      set.seed(NULL)
    } else {
      trees[[i]] = trees_i
    }
  }
  
  n.Lang <- nrow(feats)
  full.Cs <- list()
  for (n in 1:n_avg) {
    # initiate empty matrix, to be filled with languages' phylogenetic distances
    full.C.n <- matrix( data=0, nrow=n.Lang, 
                                ncol=n.Lang )
    prev <- 0
    for ( i in 1:length(trees) ){
      # trees for i^th family
      trees.i <- trees[[i]]

      # insert values of current family into full matrix
      curr <- nrow(trees.i[[1]])
      for (j in 1:curr) {
        for (k in 1:curr) {
          full.C.n[j+prev,k+prev] <- trees.i[[n]][j,k]
        }
      }
      prev <- prev + curr
    }
    full.Cs[[n]] <- full.C.n
  }
  return(full.Cs)
}

# ======== wrapper for reading both feature values and covariance matrices ========

#' suite to load data for specified families
#' @param lang character or character array, specifying language families
#' @param trait.names character array specifying trait names, at least two
#' @param avg if TRUE return average covariance matrices
#' @param n_avg number of distinct average matrices to compute if n_avg==TRUE
#' @param bin boolean, True if traits should be read as 0/1 instead of 1/2/3
#' @param other numeric, 0 or 1, if \code{bin=TRUE}, whether to read "other" category as 0/1
#' @return data in list: 
#'         N: total number of languages,
#'         r: number of traits,
#'         traits: trait values for each language as numeric array
#'         C: covariance matrix/matrices
#'         nMat: number of covariance matrices
#'         l.names: character array of names of languages 
#'                  in same order as in traits and C
#' 
get.data <- function(lang, trait.names, avg=FALSE, n_avg=1,
                     bin=FALSE, other=0){
  
  r <- length(trait.names)
  
  traits <- c()
  
  if ( length(lang)==1 ){
    feats <- read.feats( lang, bin=bin, other=other )
    for ( trait.name in trait.names ) {
      trait <- as.numeric( unlist(feats[trait.name]) )
      traits <- c(traits, trait)
    }
    lang.names <- as.character( feats$X )
  }
  
  trees.C <- read.trees(lang, lang.names, avg=avg, n_avg=n_avg)
  nMat <- length(trees.C)
  N <- nrow(trees.C[[1]])
  return(list( N=N, r=r, traits=traits, C=trees.C, nMat=nMat, l.names=lang.names ))
}

#' suite to load data for all languages
#' @param trait.names character array specifying trait names, at least two
#' @param avg if TRUE return average covariance matrices
#' @param n_avg number of distinct average matrices to compute if n_avg==TRUE
#' @param getCs if FALSE don't return covariance matrices
#' @return data in list: 
#'         Ns: array of language sizes for each family,
#'         traits: trait values for each language as numeric array ordered by family
#'         Cs: covariance matrix/matrices,
#'         nMat: number of covariance matrices
#'         l.names: character array of names of languages 
#'                  in same order as in traits and C
#' 
get.data.full <- function( trait.names, other=0, avg=TRUE, n_avg=1, getCs=TRUE){
  
  max.C <- ifelse(avg, 1, max.C)
  
  feats <- read.feats("all", bin=TRUE, other=other)
  fam.names <- unique(feats$glot)
  
  lang.names <- c()
  traits <- c()
  Ns <- c()
  for (name in fam.names) {
    fam.feats <- feats[feats$glot==name,]
    lang.names <- c( lang.names, fam.feats$X )
    Ns <- c(Ns, nrow(fam.feats))
    traits <- c( traits, as.numeric(unlist(fam.feats[trait.names[1]])) )
    traits <- c( traits, as.numeric(unlist(fam.feats[trait.names[2]])) )
  }
  
  if (getCs){
    Cs <- read.trees.multi.Lang(feats, avg=avg, n_avg=n_avg )
  } else { Cs <- NA }
  
  return(list(Ns=Ns, traits=traits, Cs=Cs, nMat=n_avg, lang.names=lang.names))
}


# ======== prepare data for cross validation ========
# eventually not used

get.data.2fold <- function( lang.name, trait.name.1, trait.name.2, 
                            M=10, p.test=3, seed.M=1, seed.N=1, fullEval=FALSE ){
  # - setup stuff: 
  `%!in%` = Negate( `%in%` )
  
  set.seed( seed.M )
  sample.M <- sample.int( M )
  set.seed( NULL )
  
  # - read traits:
  traits <- read.feats( lang.name=lang.name, bin=TRUE )
  traits$X <- as.character( traits$X )
  
  # split into folds
  N <- length( traits$X )
  N.test <- ceiling( N/p.test )
  # make sure at leas two languages are included in test data, 
  # so trees/matrices can still be read as such
  if (!fullEval){
    N.test <- max(N.test, 2)
  }
    
  
  # train-test 1
  set.seed( seed.N )
  names.test.1 <- sample( traits$X, N.test )
  indices.test.1 <- sample.int( N.test )
  set.seed( NULL )
  
  traits.test.1 <- traits[traits$X %in% names.test.1,]
  traits.train.1 <- traits[traits$X %!in% names.test.1,] 
  
  # train-test 2
  set.seed( seed.N )
  names.test.2 <- sample( traits.train.1$X, N.test )
  indices.test.2 <- sample.int( N.test )
  set.seed( NULL )
  
  traits.test.2 <- traits[traits$X %in% names.test.2,] 
  traits.train.2 <- traits[traits$X %!in% names.test.2,] 
  
  # - read trees:
  # train-test 1
  C.train.1 <- read.trees( lang.name=lang.name, feats.lang.names=traits.train.1$X )[sample.M]
  if (fullEval){
    C.test.1 <- read.trees(lang.name=lang, feats.lang.names=traits$X)[sample.M]
  } else {
    C.test.1 <- read.trees(lang.name=lang, feats.lang.names=traits.test.1$X)[sample.M]
  }
  
  # train-test 2
  C.train.2 <- read.trees( lang.name=lang.name, feats.lang.names=traits.train.2$X )[sample.M]
  if (fullEval){
    C.test.2 <- read.trees(lang.name=lang, feats.lang.names=traits$X)[sample.M]
  } else {
    C.test.2 <- read.trees(lang.name=lang, feats.lang.names=traits.test.2$X)[sample.M]
  }
  
  # - select specific traits to return
  train.1.vals <- c(traits.train.1[,trait.name.1], traits.train.1[,trait.name.2])
  if (fullEval){
    test.1.vals <- c(traits[,trait.name.1], traits[,trait.name.2])
  } else {
    test.1.vals <- c(traits.test.1[,trait.name.1], traits.test.1[,trait.name.2])
  }
  train.2.vals <- c(traits.train.2[,trait.name.1], traits.train.2[,trait.name.2])
  if (fullEval){
    test.2.vals <- c(traits[,trait.name.1], traits[,trait.name.2])
  } else {
    test.2.vals <- c(traits.test.2[,trait.name.1], traits.test.2[,trait.name.2])
  }
  test.N <- ifelse(fullEval, N, N.test)
  unseen.1 <- rep(0, N)
  unseen.1[indices.test.1] <- 1
  unseen.1 <- rep(unseen.1,2)
  unseen.2 <- rep(0, N)
  unseen.2[indices.test.2] <- 1
  unseen.2 <- rep(unseen.2,2)
  return(list(
    fold.1 = list(
      train.traits = train.1.vals,
      test.traits = test.1.vals,
      unseen = unseen.1,
      train.C = C.train.1,
      test.C = C.test.1
    ),
    fold.2 = list(
      train.traits = train.2.vals,
      test.traits = test.2.vals,
      unseen = unseen.2,
      train.C = C.train.2,
      test.C = C.test.2
    ),
    train.N = (N-N.test),
    test.N = test.N
  ))
}


get.data.3fold <- function( lang.name, trait.name.1, trait.name.2, 
                            M=10, p.test=3, seed.M=1, seed.N=1, fullEval=FALSE ){
  # - setup stuff: 
  `%!in%` = Negate( `%in%` )
  
  set.seed( seed.M )
  sample.M <- sample.int( M )
  set.seed( NULL )
  
  # - read traits:
  traits <- read.feats( lang.name=lang.name, bin=TRUE )
  traits$X <- as.character( traits$X )
  
  # split into folds
  N <- length( traits$X )
  N.test <- floor( N/p.test )
  # make sure at leas two languages are included in test data, 
  # so trees/matrices can still be read as such
  if (!fullEval){
    N.test <- max(N.test, 2)
  }
  
  
  # train-test 1
  set.seed( seed.N )
  names.test.1 <- sample( traits$X, N.test )
  indices.test.1 <- sample.int( N.test )
  set.seed( NULL )
  
  traits.test.1 <- traits[traits$X %in% names.test.1,]
  traits.train.1 <- traits[traits$X %!in% names.test.1,] 
  
  # train-test 2
  set.seed( seed.N )
  names.test.2 <- sample( traits.train.1$X, N.test )
  indices.test.2 <- sample.int( N.test )
  set.seed( NULL )
  
  traits.test.2 <- traits[traits$X %in% names.test.2,] 
  traits.train.2 <- traits[traits$X %!in% names.test.2,] 
  
  # train-test 3
  set.seed( seed.N )
  names.test.3 <- sample( traits.train.1$X[traits.train.1$X %!in% names.test.2], N.test )
  indices.test.3 <- sample.int( N.test )
  set.seed( NULL )
  
  traits.test.3 <- traits[traits$X %in% names.test.3,] 
  traits.train.3 <- traits[traits$X %!in% names.test.3,] 
  
  # - read trees:
  # train-test 1
  C.train.1 <- read.trees( lang.name=lang.name, feats.lang.names=traits.train.1$X )[sample.M]
  if (fullEval){
    C.test.1 <- read.trees(lang.name=lang, feats.lang.names=traits$X)[sample.M]
  } else {
    C.test.1 <- read.trees(lang.name=lang, feats.lang.names=traits.test.1$X)[sample.M]
  }
  
  # train-test 2
  C.train.2 <- read.trees( lang.name=lang.name, feats.lang.names=traits.train.2$X )[sample.M]
  if (fullEval){
    C.test.2 <- read.trees(lang.name=lang, feats.lang.names=traits$X)[sample.M]
  } else {
    C.test.2 <- read.trees(lang.name=lang, feats.lang.names=traits.test.2$X)[sample.M]
  }
  
  # train-test 3
  C.train.3 <- read.trees( lang.name=lang.name, feats.lang.names=traits.train.3$X )[sample.M]
  if (fullEval){
    C.test.3 <- read.trees(lang.name=lang, feats.lang.names=traits$X)[sample.M]
  } else {
    C.test.3 <- read.trees(lang.name=lang, feats.lang.names=traits.test.3$X)[sample.M]
  }
  
  
  # - select specific traits to return
  train.1.vals <- c(traits.train.1[,trait.name.1], traits.train.1[,trait.name.2])
  if (fullEval){
    test.1.vals <- c(traits[,trait.name.1], traits[,trait.name.2])
  } else {
    test.1.vals <- c(traits.test.1[,trait.name.1], traits.test.1[,trait.name.2])
  }
  train.2.vals <- c(traits.train.2[,trait.name.1], traits.train.2[,trait.name.2])
  if (fullEval){
    test.2.vals <- c(traits[,trait.name.1], traits[,trait.name.2])
  } else {
    test.2.vals <- c(traits.test.2[,trait.name.1], traits.test.2[,trait.name.2])
  }
  train.3.vals <- c(traits.train.3[,trait.name.1], traits.train.3[,trait.name.2])
  if (fullEval){
    test.3.vals <- c(traits[,trait.name.1], traits[,trait.name.2])
  } else {
    test.3.vals <- c(traits.test.3[,trait.name.1], traits.test.3[,trait.name.2])
  }
  
  test.N <- ifelse(fullEval, N, N.test)
  unseen.1 <- rep(0, N)
  unseen.1[indices.test.1] <- 1
  unseen.1 <- rep(unseen.1,2)
  unseen.2 <- rep(0, N)
  unseen.2[indices.test.2] <- 1
  unseen.2 <- rep(unseen.2,2)
  unseen.3 <- rep(0, N)
  unseen.3[indices.test.3] <- 1
  unseen.3 <- rep(unseen.3,2)
  
  return(list(
    fold.1 = list(
      train.traits = train.1.vals,
      test.traits = test.1.vals,
      unseen = unseen.1,
      train.C = C.train.1,
      test.C = C.test.1
    ),
    fold.2 = list(
      train.traits = train.2.vals,
      test.traits = test.2.vals,
      unseen = unseen.2,
      train.C = C.train.2,
      test.C = C.test.2
    ),
    fold.3 = list(
      train.traits = train.3.vals,
      test.traits = test.3.vals,
      unseen = unseen.3,
      train.C = C.train.3,
      test.C = C.test.3
    ),
    train.N = (N-N.test),
    test.N = test.N
  ))
}


