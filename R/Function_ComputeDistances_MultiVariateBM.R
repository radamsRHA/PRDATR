#' Function_ComputeDistances_MultiVariateBM: function to compute the Hellinger and Kl distance for two probabilisty phylogenetic models under continuous trait evolution
#'
#' This function returns a vector containing the Hellinger and KL distances between two tree models
#' @param list.Model_01 List containing the following (1) handle.Phylogeny, (2) vector.Z = vector of mean (ancestral) state values, and (3) vector.Model_01_Theta = vector containing relevant parameters for the models, and (4) matrix.C = matrix containing the across trait variance/covariances
#' @param list.Model_02 List containing the following (1) handle.Phylogeny, (3) vector.Z = vector of mean (ancestral) state values, and (3) vector.Model_02_Theta = vector containing relevant parameters for the models, and (4) matrix.C = matrix containing the across trait variance/covariances
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
#' #################################################
#' # Set a vector containing parameters for Models #
#' #################################################
#' vector.Model_01_Theta <- c(1)
#' names(vector.Model_01_Theta) <- c("Sig2")
#' vector.Model_02_Theta <- c(1)
#' names(vector.Model_02_Theta) <- c("Sig2")
#'
#' matrix.Model_01_C <- matrix(nrow = 2, ncol = 2)
#' matrix.Model_01_C[1,] <- c(1, 0.5)
#' matrix.Model_01_C[2,] <- c(0.5, 1)
#'
#' matrix.Model_02_C <- matrix.Model_01_C
#' matrix.Model_02_C[1,2] <- matrix.Model_02_C[2,1] <- 0.75
#'
#' #####################
#' # Build Model lists #
#' #####################
#'
#' list.Model_01 <- list(handle.Phylogeny = handle.Figure01_Felsenstein1985_Tree,
#'                       vector.Z = rep(0, 2*length(handle.Figure01_Felsenstein1985_Tree$tip.label)),
#'                       vector.Theta = vector.Model_01_Theta,
#'                       matrix.C =  matrix.Model_01_C)
#'
#' list.Model_02 <- list(handle.Phylogeny = handle.Figure01_Felsenstein1985_Tree,
#'                       vector.Z = rep(0, 2*length(handle.Figure01_Felsenstein1985_Tree$tip.label)),
#'                       vector.Theta = vector.Model_02_Theta,
#'                       matrix.C =  matrix.Model_02_C)
#'
#' #####################
#' # Compute distances #
#' #####################
#' Function_ComputeDistances_MultiVariateBM(list.Model_01 = list.Model_01, list.Model_02 = list.Model_02)
#'
#'



############################################
# Function_ComputeDistances_MultiVariateBM #
############################################
Function_ComputeDistances_MultiVariateBM <- function(list.Model_01, list.Model_02){

  ##################
  # Build Model 01 #
  ##################
  handle.Model_01_Tree <- list.Model_01$handle.Phylogeny
  vector.Model_01_Z <- list.Model_01$vector.Z
  vector.Model_01_Theta <- list.Model_01$vector.Theta
  numeric.Model_01_Sig2 <- vector.Model_01_Theta["Sig2"]
  handle.Model_01_ReScaled <- matrix.Model_01_VarCovar <- NULL
  matrix.Model_01_C <- list.Model_01$matrix.C

  ########################
  # Get vcv for Model 01 #
  ########################
  handle.Model_01_ReScaled <- rescale(x = handle.Model_01_Tree, model = "BM")
  handle.Model_01_ReScaled <- handle.Model_01_ReScaled(sigsq = numeric.Model_01_Sig2)
  matrix.Model_01_VarCovar <- vcv(phy = handle.Model_01_ReScaled)

  #####################################
  # Get multivariate vcv for Model 01 #
  #####################################
  #matrix.Model_01_VarCovar_MultiVariate <- kronecker(X = matrix.Model_01_VarCovar, Y = matrix.Model_01_C)
  matrix.Model_01_VarCovar_MultiVariate <- kronecker(X = matrix.Model_01_C, Y = matrix.Model_01_VarCovar)

  ##################
  # Build Model 02 #
  ##################
  handle.Model_02_Tree <- list.Model_02$handle.Phylogeny
  vector.Model_02_Z <- list.Model_02$vector.Z
  vector.Model_02_Theta <- list.Model_02$vector.Theta
  numeric.Model_02_Sig2 <- vector.Model_02_Theta["Sig2"]
  handle.Model_02_ReScaled <- matrix.Model_02_VarCovar <- NULL
  matrix.Model_02_C <- list.Model_02$matrix.C

  ###########################################
  # Specify evolutionary model for Model 02 #
  ###########################################
  handle.Model_02_ReScaled <- rescale(x = handle.Model_02_Tree, model = "BM")
  handle.Model_02_ReScaled <- handle.Model_02_ReScaled(sigsq = numeric.Model_02_Sig2)
  matrix.Model_02_VarCovar <- vcv(phy = handle.Model_02_ReScaled)

  #####################################
  # Get multivariate vcv for Model 02 #
  #####################################
  #matrix.Model_02_VarCovar_MultiVariate <- kronecker(X = matrix.Model_02_VarCovar, Y = matrix.Model_02_C)
  matrix.Model_02_VarCovar_MultiVariate <- kronecker(X = matrix.Model_02_C, Y = matrix.Model_02_VarCovar)


  #################################
  # Compute Hellinger coefficient #
  #################################
  numeric.HellingerCoefficient <- normdiff(mu1=vector.Model_01_Z,sigma1=matrix.Model_01_VarCovar_MultiVariate,
                                           mu2=vector.Model_02_Z,sigma2=matrix.Model_02_VarCovar_MultiVariate,
                                           method="Hellinger", inv = F, s = 0.5)[1]

  numeric.Distance_Hellinger <- 1-numeric.HellingerCoefficient

  ###################################
  # Compute KL distance coefficient #
  ###################################
  numeric.Distance_KL <- normdiff(mu1=vector.Model_01_Z,sigma1=matrix.Model_01_VarCovar_MultiVariate,
                                  mu2=vector.Model_02_Z,sigma2=matrix.Model_02_VarCovar_MultiVariate,
                                  method="KL")[1]

  vector.Distances <- c(numeric.Distance_Hellinger, numeric.Distance_KL)
  names(vector.Distances) <- c("dH", "dKL")

  return(vector.Distances)
}


