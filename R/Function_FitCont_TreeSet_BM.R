#' Function_FitCont_TreeSet_BM: function to estimate BM model parameters for each tree in a given set
#'
#' This function returns a matrix of fitted model parameters (assuming a BM) for a given set of trees
#' @param phylogeny.TreeSet List of phylogenetic trees
#' @param vector.InputData Vector of continuous trait data
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
#'
#' #########################
#' # Simulate random trees #
#' #########################
#' handle.SimulatedRandomTrees <- rmtree(N = 10, n = 10)
#' vector.SimulatedDataset <- fastBM(tree = handle.SimulatedRandomTrees[[1]], a = 0, mu = 1)
#'
#' handle.FittedPhylogeneticModels <- Function_FitCont_TreeSet_BM(phylogeny.TreeSet = handle.SimulatedRandomTrees, vector.InputData = vector.SimulatedDataset)
#'
###############################
# Function_FitCont_TreeSet_BM #
###############################
Function_FitCont_TreeSet_BM <- function(phylogeny.TreeSet, vector.InputData){

  ###################
  # Summarize input #
  ###################
  numeric.NumberOfTrees <- length(phylogeny.TreeSet)
  matrix.ParamEstimates <- matrix(nrow = numeric.NumberOfTrees, ncol = 3)
  colnames(matrix.ParamEstimates) <- c("z", "sig2", "lnl")

  ######################
  # Loop through trees #
  ######################
  for (i in 1:numeric.NumberOfTrees){
    print(i)

    ################
    # Extract Tree #
    ################
    phylogeny.Tree_i <- phylogeny.TreeSet[[i]]
    handle.Results_FitCont <- fitContinuous(phy = phylogeny.Tree_i, dat = vector.InputData, SE = 0, model = "BM")
    matrix.ParamEstimates[i, "z"] <- handle.Results_FitCont$opt$z0
    matrix.ParamEstimates[i, "sig2"] <- handle.Results_FitCont$opt$sigsq
    matrix.ParamEstimates[i, "lnl"] <- handle.Results_FitCont$opt$lnL


  }

  return(matrix.ParamEstimates)
}
