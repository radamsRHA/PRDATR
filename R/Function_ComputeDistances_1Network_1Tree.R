#' Function_ComputeDistances_1Network_1Tree: function to compute probabilistic distances between a network and a bifurcating tree
#'
#' This function returns a vector containing the Hellinger and KL distances between two tree models
#' @param list.Model_01_Network List containing the following (1) handle.Phylogeny, (2) string.Model = "BM", "OU", "EB", "lambda", "kappa", "delta", (3) vector.Z = vector of mean (ancestral) state values, and (4) vector.Model_01_Theta = vector containing relevant parameters for the models
#' @param list.Model_02_Tree List containing the following (1) handle.Phylogeny, (2) string.Model = "BM", "OU", "EB", "lambda", "kappa", "delta", (3) vector.Z = vector of mean (ancestral) state values, and (4) vector.Model_01_Theta = vector containing relevant parameters for the models
#' @param boo.SortNames Boolean (TRUE or FALSE; default = FALSE). Will order the VCV for the second model based on the rownames and colnames of the VCV generated for the first model. Try this option if using different trees for the two models
#' @keywords probabilistic phylogenetic distances, model distance, brownian motion, continuous trait evolution
#' @return vector.Distances Vector containing the distances computed between the two focal tree models\cr
#' @export
#' @examples
#'
#'
#' ################
#' # Load depends #
#' ################
#' library(ape)
#' library(geiger)
#' library(gaussDiff)
#' library(BMhyd)
#'
#' ######################
#' # Simulate a network #
#' ######################
#' handle.SimulatedNetwork <- SimulateNetwork(ntax.nonhybrid=5, ntax.hybrid=3,
#'                                            flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
#'                                            tree.height = 1, allow.ghost=FALSE)
#'
#' ###############################
#' # Set the first network model #
#' ###############################
#' list.Model01_BM <- list(handle.Phylogeny = handle.SimulatedNetwork$phy,
#'                         handle.Flow = handle.SimulatedNetwork$flow,
#'                         string.Model = "BM",
#'                         vector.Z = rep(0, 8),
#'                         numeric.Sig2 = 1)
#'
#' #############################
#' # Set the second tree model #
#' #############################
#' list.Model02_BM <- list(handle.Phylogeny = handle.SimulatedNetwork$phy,
#'                         string.Model = "BM",
#'                         vector.Z = rep(0, 8),
#'                         numeric.Sig2 = 1)
#'
#'
#' #####################
#' # Compute distances #
#' #####################
#' Function_ComputeDistances_1Network_1Tree(list.Model_01_Network = list.Model01_BM, list.Model_02_Tree = list.Model02_BM)
#'

#############################################
# Function_ComputeDistances_1Network_1Tree #
#############################################
Function_ComputeDistances_1Network_1Tree <- function(list.Model_01_Network, list.Model_02_Tree, boo.SortNames = T){


  ##################
  # Build Model 01 #
  ##################
  handle.Model_01_Tree <- list.Model_01_Network$handle.Phylogeny
  handle.Model_01_Flow <- list.Model_01_Network$handle.Flow
  vector.Model_01_Z <- list.Model_01_Network$vector.Z
  numeric.Model_01_Sig2 <- list.Model_01_Network$numeric.Sig2
  handle.Model_01_ReScaled <- matrix.Model_01_VarCovar <- NULL
  matrix.VFV_Model_01_Network <- GetVModified(x = c(numeric.Model_01_Sig2,vector.Model_01_Z[1],SE = 0),
                                              phy = handle.Model_01_Tree,
                                              flow = handle.Model_01_Flow, "bt")

  ##################
  # Build Model 02 #
  ##################
  handle.Model_02_Tree <- list.Model_02_Tree$handle.Phylogeny
  vector.Model_02_Z <- list.Model_02_Tree$vector.Z
  numeric.Model_02_Sig2 <- list.Model_02_Tree$numeric.Sig2
  handle.Model_02_ReScaled <- matrix.Model_02_VarCovar <- NULL

  handle.Model_02_Tree_ReScaled <- rescale(x = handle.Model_02_Tree, model = "BM")
  handle.Model_02_Tree_ReScaled <- handle.Model_02_Tree_ReScaled(sigsq = numeric.Model_02_Sig2)
  matrix.Model_02_Tree_VarCovar <- vcv(phy = handle.Model_02_Tree_ReScaled)
  
  #################################  
  # Sort matrices by the firs one #
  #################################
  if (boo.SortNames == TRUE){
    matrix.Model_02_Tree_VarCovar <- matrix.Model_02_Tree_VarCovar[rownames(matrix.VFV_Model_01_Network), colnames(matrix.VFV_Model_01_Network)]
   }

  ##############################
  # Compute Hellinger distance #
  ##############################
  numeric.HellingerCoefficient <- normdiff(mu1=vector.Model_01_Z, sigma1= matrix.VFV_Model_01_Network,
                                           mu2=vector.Model_02_Z, sigma2= matrix.Model_02_Tree_VarCovar,
                                           method="Hellinger", inv = F, s = 0.5)[1]

  numeric.Distance_Hellinger <- 1-numeric.HellingerCoefficient


  ###################################
  # Compute KL distance coefficient #
  ###################################
  numeric.Distance_KL <- normdiff(mu1=vector.Model_01_Z, sigma1=matrix.VFV_Model_01_Network,
                                  mu2=vector.Model_02_Z, sigma2=matrix.Model_02_Tree_VarCovar,
                                  method="KL")[1]

  vector.Distances <- c(numeric.Distance_Hellinger, numeric.Distance_KL)
  names(vector.Distances) <- c("dH", "dKL")

  return(vector.Distances)
}
