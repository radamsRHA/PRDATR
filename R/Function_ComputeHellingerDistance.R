#' Function_ComputeHellingerDistances: function to compute the Hellinger and Kl distance for two probabilisty phylogenetic models under continuous trait evolution
#'
#' This function returns a vector containing the Hellinger and KL distances between two tree models
#' @param list.Model_01 List containing the following (1) handle.Phylogeny, (2) string.Model = "BM", "OU", "EB", "lambda", "kappa", "delta", (3) vector.Z = vector of mean (ancestral) state values, and (4) vector.Model_01_Theta = vector containing relevant parameters for the models
#' @param list.Model_02 List containing the following (1) handle.Phylogeny, (2) string.Model = "BM", "OU", "EB", "lambda", "kappa", "delta", (3) vector.Z = vector of mean (ancestral) state values, and (4) vector.Model_01_Theta = vector containing relevant parameters for the models
#' @param boo.SortNames Boolean (TRUE or FALSE; default = FALSE). Will order the VCV for the second model based on the rownames and colnames of the VCV generated for the first model. Try this option if using different trees for the two models
#' @keywords probabilistic phylogenetic distances, model distance, brownian motion, continuous trait evolution
#' @return vector.Distances Vector containing the distances computed between the two focal tree models\cr
#' @export
#' @examples
#'
#' ################
#' # Load depends #
#' ################
#' library(ape)
#' library(geiger)
#' library(gaussDiff)
#'
#' ############################################################################################
#' # Specifity example tree (Fig. 8, Felsenstein 1985) used for demonstrating model distances #
#' ############################################################################################
#' string.Figure01_Felsenstein1985_Tree <- "(((Species_7:1.635983031,Species_8:0.8079384444):1.801052391,(Species_6:0.4510394335,Species_5:1.208543249):0.4146682462):0.4476024657,((Species_4:0.9434278477,Species_3:0.2480806489):2.642993143,(Species_2:2.588686978,Species_1:0.908678702):0.3223280154):1.399150111):1;"
#' handle.Figure01_Felsenstein1985_Tree <- read.tree(text = string.Figure01_Felsenstein1985_Tree)
#'
#' ###################################################
#' # Set a vector containing parameters for Model 01 #
#' ###################################################
#' vector.Model_02_Theta <- c(1, 1)
#' names(vector.Model_02_Theta) <- c("Sig2", "alpha")
#' vector.Model_01_Theta <- c(1)
#' names(vector.Model_01_Theta) <- c("Sig2")
#'
#'
#' list.Model01_BM <- list(handle.Phylogeny = handle.Figure01_Felsenstein1985_Tree,
#'                         string.Model = "BM",
#'                         vector.Z = rep(0, length(handle.Figure01_Felsenstein1985_Tree$tip.label)),
#'                         vector.Theta = vector.Model_01_Theta)
#'
#' list.Model02_BM <- list(handle.Phylogeny = handle.Figure01_Felsenstein1985_Tree,
#'                         string.Model = "OU",
#'                         vector.Z = rep(0, length(handle.Figure01_Felsenstein1985_Tree$tip.label)),
#'                         vector.Theta = vector.Model_02_Theta)
#'
#' Function_ComputeHellingerDistances(list.Model_01 = list.Model01_BM, list.Model_02 = list.Model02_BM)
#'

######################################
# Function_ComputeHellingerDistances #
######################################
Function_ComputeHellingerDistances <- function(list.Model_01, list.Model_02, boo.SortNames = NA){

  ##################
  # Build Model 01 #
  ##################
  handle.Model_01_Tree <- list.Model_01$handle.Phylogeny
  string.Model_01_Model <- list.Model_01$string.Model
  vector.Model_01_Z <- list.Model_01$vector.Z
  vector.Model_01_Theta <- list.Model_01$vector.Theta
  numeric.Model_01_Sig2 <- vector.Model_01_Theta["Sig2"]
  handle.Model_01_ReScaled <- matrix.Model_01_VarCovar <- NULL

  ###########################################
  # Specify evolutionary model for Model 01 #
  ###########################################
  if (string.Model_01_Model == "BM"){
    handle.Model_01_ReScaled <- rescale(x = handle.Model_01_Tree, model = "BM")
    handle.Model_01_ReScaled <- handle.Model_01_ReScaled(sigsq = numeric.Model_01_Sig2)
    matrix.Model_01_VarCovar <- vcv(phy = handle.Model_01_ReScaled)
  }
  if (string.Model_01_Model == "OU"){
    handle.Model_01_ReScaled <- rescale(x = handle.Model_01_Tree, model = "OU")
    handle.Model_01_ReScaled <- handle.Model_01_ReScaled(alpha = vector.Model_01_Theta["alpha"], sigsq = numeric.Model_01_Sig2)
    matrix.Model_01_VarCovar <- vcv(phy = handle.Model_01_ReScaled)
  }
  if (string.Model_01_Model == "EB"){
    handle.Model_01_ReScaled <- rescale(x = handle.Model_01_Tree, model = "EB")
    handle.Model_01_ReScaled <- handle.Model_01_ReScaled(a = vector.Model_01_Theta["a"], sigsq = numeric.Model_01_Sig2)
    matrix.Model_01_VarCovar <- vcv(phy = handle.Model_01_ReScaled)
  }
  if (string.Model_01_Model == "lambda"){
    handle.Model_01_ReScaled <- rescale(x = handle.Model_01_Tree, model = "lambda")
    handle.Model_01_ReScaled <- handle.Model_01_ReScaled(lambda = vector.Model_01_Theta["lambda"], sigsq = numeric.Model_01_Sig2)
    matrix.Model_01_VarCovar <- vcv(phy = handle.Model_01_ReScaled)
  }
  if (string.Model_01_Model == "kappa"){
    handle.Model_01_ReScaled <- rescale(x = handle.Model_01_Tree, model = "kappa")
    handle.Model_01_ReScaled <- handle.Model_01_ReScaled(kappa = vector.Model_01_Theta["kappa"], sigsq = numeric.Model_01_Sig2)
    matrix.Model_01_VarCovar <- vcv(phy = handle.Model_01_ReScaled)
  }
  if (string.Model_01_Model == "delta"){
    handle.Model_01_ReScaled <- rescale(x = handle.Model_01_Tree, model = "delta")
    handle.Model_01_ReScaled <- handle.Model_01_ReScaled(delta = vector.Model_01_Theta["delta"], sigsq = numeric.Model_01_Sig2)
    matrix.Model_01_VarCovar <- vcv(phy = handle.Model_01_ReScaled)
  }
  if (string.Model_01_Model == "ACDC"){
    #handle.Model_01_Tree$edge.length <- handle.Model_01_Tree$edge.length/max(handle.Model_01_Tree$edge.length)
    matrix.Model_01_VarCovar <- vcv(phy = handle.Model_01_Tree)

    numeric.r <- vector.Model_01_Theta["r"]
    numeric.Sig2 <- numeric.Model_01_Sig2
    numeric.alpha <- numeric.r/2
    numeric.Sig2.0 <- numeric.Sig2*exp(-2*numeric.alpha*1)

    #matrix.Model_01_VarCovar <- round(x = numeric.Sig2.0*(exp(numeric.r*matrix.Model_01_VarCovar)-1)/numeric.r, digits = 4)
    matrix.Model_01_VarCovar <- numeric.Sig2.0*(exp(numeric.r*matrix.Model_01_VarCovar)-1)/numeric.r
  }

  ##################
  # Build Model 02 #
  ##################
  handle.Model_02_Tree <- list.Model_02$handle.Phylogeny
  string.Model_02_Model <- list.Model_02$string.Model
  vector.Model_02_Z <- list.Model_02$vector.Z
  vector.Model_02_Theta <- list.Model_02$vector.Theta
  numeric.Model_02_Sig2 <- vector.Model_02_Theta["Sig2"]
  handle.Model_02_ReScaled <- matrix.Model_02_VarCovar <- NULL

  ###########################################
  # Specify evolutionary model for Model 02 #
  ###########################################
  if (string.Model_02_Model == "BM"){
    handle.Model_02_ReScaled <- rescale(x = handle.Model_02_Tree, model = "BM")
    handle.Model_02_ReScaled <- handle.Model_02_ReScaled(sigsq = numeric.Model_02_Sig2)
    matrix.Model_02_VarCovar <- vcv(phy = handle.Model_02_ReScaled)
  }
  if (string.Model_02_Model == "OU"){
    handle.Model_02_ReScaled <- rescale(x = handle.Model_02_Tree, model = "OU")
    handle.Model_02_ReScaled <- handle.Model_02_ReScaled(alpha = vector.Model_02_Theta["alpha"], sigsq = numeric.Model_02_Sig2)
    matrix.Model_02_VarCovar <- vcv(phy = handle.Model_02_ReScaled)
  }
  if (string.Model_02_Model == "EB"){
    handle.Model_02_ReScaled <- rescale(x = handle.Model_02_Tree, model = "EB")
    handle.Model_02_ReScaled <- handle.Model_02_ReScaled(a = vector.Model_02_Theta["a"], sigsq = numeric.Model_02_Sig2)
    matrix.Model_02_VarCovar <- vcv(phy = handle.Model_02_ReScaled)
  }
  if (string.Model_02_Model == "lambda"){
    handle.Model_02_ReScaled <- rescale(x = handle.Model_02_Tree, model = "lambda")
    handle.Model_02_ReScaled <- handle.Model_02_ReScaled(lambda = vector.Model_02_Theta["lambda"], sigsq = numeric.Model_02_Sig2)
    matrix.Model_02_VarCovar <- vcv(phy = handle.Model_02_ReScaled)
  }
  if (string.Model_02_Model == "kappa"){
    handle.Model_02_ReScaled <- rescale(x = handle.Model_02_Tree, model = "kappa")
    handle.Model_02_ReScaled <- handle.Model_02_ReScaled(kappa = vector.Model_02_Theta["kappa"], sigsq = numeric.Model_02_Sig2)
    matrix.Model_02_VarCovar <- vcv(phy = handle.Model_02_ReScaled)
  }
  if (string.Model_02_Model == "delta"){
    handle.Model_02_ReScaled <- rescale(x = handle.Model_02_Tree, model = "delta")
    handle.Model_02_ReScaled <- handle.Model_02_ReScaled(delta = vector.Model_02_Theta["delta"], sigsq = numeric.Model_02_Sig2)
    matrix.Model_02_VarCovar <- vcv(phy = handle.Model_02_ReScaled)
  }
  if (string.Model_02_Model == "ACDC"){
    #handle.Model_02_Tree$edge.length <- handle.Model_02_Tree$edge.length/max(handle.Model_02_Tree$edge.length)
    matrix.Model_02_VarCovar <- vcv(phy = handle.Model_02_Tree)

    numeric.r <- vector.Model_02_Theta["r"]
    numeric.Sig2 <- numeric.Model_02_Sig2
    numeric.alpha <- numeric.r/2
    numeric.Sig2.0 <- numeric.Sig2*exp(-2*numeric.alpha*1)

    #matrix.Model_02_VarCovar <- round(x = numeric.Sig2.0*(exp(numeric.r*matrix.Model_02_VarCovar)-1)/numeric.r, digits = 4)
    matrix.Model_02_VarCovar <- numeric.Sig2.0*(exp(numeric.r*matrix.Model_02_VarCovar)-1)/numeric.r
  }
  #################################  
  # Sort matrices by the firs one #
  #################################
  if (boo.SortNames == TRUE){
    matrix.Model_02_VarCovar <- matrix.Model_02_VarCovar[rownames(matrix.Model_01_VarCovar), colnames(matrix.Model_01_VarCovar)]
   }

  #################################
  # Compute Hellinger coefficient #
  #################################
  numeric.HellingerCoefficient <- normdiff(mu1=vector.Model_01_Z,sigma1=matrix.Model_01_VarCovar,
                                           mu2=vector.Model_02_Z,sigma2=matrix.Model_02_VarCovar,
                                           method="Hellinger", inv = F, s = 0.5)[1]

  numeric.Distance_Hellinger <- 1-numeric.HellingerCoefficient
  names(numeric.Distance_Hellinger) <- "dH"

  return(numeric.Distance_Hellinger)
}


