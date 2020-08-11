#' Function_ComputePairWiseModelDistances_BM: function to compute pairwise model distances between a set of trees given their model parameter estimates
#'
#' This function returns a list of pairwise model distances
#' @param phylogeny.TreeSet Set of phylogenetic trees
#' @param matrix.ModelParams Matrix of model parameters (BM) for each tree in the given set
#' @keywords probabilistic phylogenetic distances, model distance, brownian motion, continuous trait evolution
#' @return List List of pairwise model distances
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
#' #########################
#' # Simulate random trees #
#' #########################
#' handle.SimulatedRandomTrees <- rmtree(N = 10, n = 10)
#' vector.SimulatedDataset <- fastBM(tree = handle.SimulatedRandomTrees[[1]], a = 0, mu = 1)
#'
#' #######################
#' # Fit models to trees #
#' #######################
#' handle.FittedPhylogeneticModels <- Function_FitCont_TreeSet_BM(phylogeny.TreeSet = handle.SimulatedRandomTrees, vector.InputData = vector.SimulatedDataset)
#'
#' ##############################
#' # Compute pairwise distances #
#' ##############################
#' handle.Results <- Function_ComputePairWiseModelDistances_BM(phylogeny.TreeSet = handle.SimulatedRandomTrees, matrix.ModelParams = handle.FittedPhylogeneticModels)
#'

#############################################
# Function_ComputePairWiseModelDistances_BM #
#############################################
Function_ComputePairWiseModelDistances_BM <- function(phylogeny.TreeSet, matrix.ModelParams){

  ###################
  # Summarize input #
  ###################
  numeric.NumberOfTreeModels <- length(phylogeny.TreeSet)
  matrix.PairwiseDistances_H <- matrix(nrow = numeric.NumberOfTreeModels, ncol = numeric.NumberOfTreeModels)
  diag(matrix.PairwiseDistances_H) <- 0
  matrix.PairwiseDistances_KL <- matrix(nrow = numeric.NumberOfTreeModels, ncol = numeric.NumberOfTreeModels)
  diag(matrix.PairwiseDistances_KL) <- 0

  ############################
  # Loop through tree models #
  ############################
  for (i in 1:numeric.NumberOfTreeModels){

    ##########################
    # Extract model for tree #
    ##########################
    phylogeny.Tree_01 <- phylogeny.TreeSet[[i]]
    numeric.n_01 <- length(phylogeny.Tree_01$tip.label)
    vector.ModelParams_01 <- matrix.ModelParams[i,]
    numeric.z_01 <- vector.ModelParams_01["z"]
    numeric.Sig2_01 <- vector.ModelParams_01["sig2"]
    vector.Theta_01 <- c(numeric.Sig2_01)
    names(vector.Theta_01) <- c("Sig2")

    handle.List_Model_01 <- list(handle.Phylogeny = phylogeny.Tree_01,
                                 string.Model = "BM",
                                 vector.Z = rep(numeric.z_01, numeric.n_01),
                                 vector.Theta = vector.Theta_01)


    ############################
    # Loop through tree models #
    ############################
    for (j in 1:numeric.NumberOfTreeModels){

      if (j!=i){

        ##########################
        # Extract model for tree #
        ##########################
        phylogeny.Tree_02 <- phylogeny.TreeSet[[j]]
        numeric.n_02 <- length(phylogeny.Tree_02$tip.label)
        vector.ModelParams_02 <- matrix.ModelParams[j,]
        numeric.z_02 <- vector.ModelParams_02["z"]
        numeric.Sig2_02 <- vector.ModelParams_02["sig2"]
        vector.Theta_02 <- c(numeric.Sig2_02)
        names(vector.Theta_02) <- c("Sig2")

        print(c(i,j))
        handle.List_Model_02 <- list(handle.Phylogeny = phylogeny.Tree_02,
                                     string.Model = "BM",
                                     vector.Z = rep(numeric.z_02, numeric.n_02),
                                     vector.Theta = vector.Theta_02)

        ####################
        # Compute distance #
        ####################
        vector.Distances <- Function_ComputeDistances(list.Model_01 = handle.List_Model_01,
                                                      list.Model_02 = handle.List_Model_02)

        matrix.PairwiseDistances_H[i,j] <- vector.Distances["dH"]
        matrix.PairwiseDistances_KL[i,j] <- vector.Distances["dKL"]




      }
    }
  }
  return(list(matrix.PairwiseDistances_H = matrix.PairwiseDistances_H, matrix.PairwiseDistances_KL = matrix.PairwiseDistances_KL,
              vector.PairwiseDistances_H = matrix.PairwiseDistances_H[upper.tri(matrix.PairwiseDistances_H)]))
}
