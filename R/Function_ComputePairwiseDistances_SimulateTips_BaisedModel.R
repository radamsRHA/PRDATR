#' Function_ComputePairwiseDistances_SimulateTips_BaisedModel: function to generate pairwise distance computations
#'
#' This function returns a list with (1) matrix.PairwiseDistances_H, (2) matrix.PairwiseDistances_KL, and (3)vector.PairwiseDistances_H
#' @param numeric.NumberReps Number of simulation replicates
#' @param numeric.n Number of tips
#' @param numeric.birth Birth Rate
#' @param string.Model String for evolutonary model
#' @keywords probabilistic phylogenetic distances, model distance, brownian motion, continuous trait evolution
#' @return List list with (1) matrix.PairwiseDistances_H, (2) matrix.PairwiseDistances_KL, and (3)vector.PairwiseDistances_H
#' @export
#' @examples
#'
#'
#'

##############################################################
# Function_ComputePairwiseDistances_SimulateTips_BaisedModel #
##############################################################
Function_ComputePairwiseDistances_SimulateTips_BaisedModel <- function(numeric.NumberReps, numeric.n, numeric.bias, string.Model, vector.Theta){

  ###################
  # Summarize input #
  ###################
  matrix.PairwiseDistances_H <- matrix(nrow = numeric.NumberReps, ncol = numeric.NumberReps)
  diag(matrix.PairwiseDistances_H) <- 0
  matrix.PairwiseDistances_KL <- matrix(nrow = numeric.NumberReps, ncol = numeric.NumberReps)
  diag(matrix.PairwiseDistances_KL) <- 0

  #########################################
  # Simulate a list of phylogenetic trees #
  #########################################
  list.SimulatedTrees <- list()
  class(list.SimulatedTrees) <- "multiPhylo"
  
  for (t in 1:numeric.NumberReps){
    
    list.SimulatedTrees[[t]] <- as.phylo.treeshape(x = rbiased(tip.number = numeric.n, p = numeric.bias))
    #handle.Tree <- as.phylo.treeshape(x = rbiased(tip.number = numeric.n, p = numeric.bias))
    #handle.Tree <- rescale(x = handle.Tree, model = "depth")
    #handle.Tree <- handle.Tree(depth = 1)
    #list.SimulatedTrees[[t]] <- handle.Tree
    
  }

  ################################
  # Loop through simulated trees #
  ################################
  for (i in 1:numeric.NumberReps){

    ##################
    # Extract tree 1 #
    ##################
    handle.Tree_1 <- list.SimulatedTrees[[i]]
    handle.List_Model_01 <- list(handle.Phylogeny = handle.Tree_1,
                                 string.Model = string.Model,
                                 vector.Z = rep(0, numeric.n),
                                 vector.Theta = vector.Theta)

    ################################
    # Loop through simulated trees #
    ################################
    for (j in 1:numeric.NumberReps){


      ##################
      # Extract Tree 2 #
      ##################
      handle.Tree_2 <- list.SimulatedTrees[[j]]

      ##############################
      # Compute off-diag distances #
      ##############################
      if (i !=j){


        ##################
        # Build Model 02 #
        ##################
        handle.List_Model_02 <- list(handle.Phylogeny = handle.Tree_2,
                                     string.Model = string.Model,
                                     vector.Z = rep(0, numeric.n),
                                     vector.Theta = vector.Theta)
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
