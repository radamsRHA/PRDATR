#' Function_ComputeDistances_2Networks: function to compute probabilistic distances between two network models
#'
#' This function returns a vector containing the Hellinger and KL distances between two tree models
#' @param list.Model_01_Network List containing network model with (1) phylogenetic tree (2) table of migratio flow (3) vector.Z and (4) Sig2 parameter
#' @param list.Model_02_Network List containing network model with (1) phylogenetic tree (2) table of migratio flow (3) vector.Z and (4) Sig2 parameter
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
#' ######################################
#' # Specific list for network model 01 #
#' ######################################
#' list.Model_01_Network <- list(handle.Phylogeny = handle.SimulatedNetwork$phy,
#'                               handle.Flow = handle.SimulatedNetwork$phy,
#'                               vector.Z = rep(0, length(network$phy$tip.label)),
#'                               numeric.Sig2 = 1)
#'
#' list.Model_02_Network <- list(handle.Phylogeny = handle.SimulatedNetwork$phy,
#'                               handle.Flow = handle.SimulatedNetwork$phy,
#'                               vector.Z = rep(0, length(network$phy$tip.label)),
#'                               numeric.Sig2 = 2)
#'
#' #####################
#' # Compute distances #
#' #####################
#' Function_ComputeDistances_2Networks(list.Model_01_Network = list.Model_01_Network,
#'                                     list.Model_02_Network = list.Model_02_Network)
#'

#######################################
# Function_ComputeDistances_2Networks #
#######################################
Function_ComputeDistances_2Networks <- function(list.Model_01_Network, list.Model_02_Network){

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
  handle.Model_02_Tree <- list.Model_02_Network$handle.Phylogeny
  handle.Model_02_Flow <- list.Model_02_Network$handle.Flow
  vector.Model_02_Z <- list.Model_02_Network$vector.Z
  numeric.Model_02_Sig2 <- list.Model_02_Network$numeric.Sig2
  handle.Model_02_ReScaled <- matrix.Model_02_VarCovar <- NULL
  matrix.VFV_Model_02_Network <- GetVModified(x = c(numeric.Model_02_Sig2,vector.Model_02_Z[1],SE = 0),
                                          phy = handle.Model_02_Tree,
                                          flow = handle.Model_02_Flow, "bt")

  ##############################
  # Compute Hellinger distance #
  ##############################
  numeric.HellingerCoefficient <- normdiff(mu1=vector.Model_01_Z, sigma1= matrix.VFV_Model_01_Network,
                                           mu2=vector.Model_02_Z, sigma2= matrix.VFV_Model_02_Network,
                                           method="Hellinger", inv = F, s = 0.5)[1]

  numeric.Distance_Hellinger <- 1-numeric.HellingerCoefficient


  ###################################
  # Compute KL distance coefficient #
  ###################################
  numeric.Distance_KL <- normdiff(mu1=vector.Model_01_Z, sigma1=matrix.VFV_Model_01_Network,
                                  mu2=vector.Model_02_Z, sigma2=matrix.VFV_Model_02_Network,
                                  method="KL")[1]

  vector.Distances <- c(numeric.Distance_Hellinger, numeric.Distance_KL)
  names(vector.Distances) <- c("dH", "dKL")

  return(vector.Distances)
}
