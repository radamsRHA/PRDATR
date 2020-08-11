#' Function_ConductSimExperimentOU_ModelTestBM_OU: function to compute p-values for likelihood ratio test and compute distances
#'
#' This function returns a matrix.SimulationResults for the simulation results
#' @param list.Model_Simulation_OU List containing components of simulation model
#' @keywords probabilistic phylogenetic distances, model distance, brownian motion, continuous trait evolution
#' @return vector.Distances Vector containing the distances computed between the two focal tree models\cr
#' @export
#' @examples
#'
#' ################
#' # Load depends #
#' ################
#' #library(ape)
#' #library(phytools)
#' #library(geiger)
#' #library(vioplot)
#' #library(gaussDiff)
#' #library(MASS)
#'
#' ##########################
#' # Build simulation model #
#' ##########################
#' phylogeny.RandomTree <- pbtree(b = 1, d = 0, n = 50)
#' vector.Model_Simulation_Theta <- c(1, 10)
#' names(vector.Model_Simulation_Theta) <- c("Sig2", "alpha")
#'
#' list.Model_Simulation_OU <- list(handle.Phylogeny = phylogeny.RandomTree,
#'                                  string.Model = "OU",
#'                                  vector.Z = rep(0, length(phylogeny.RandomTree$tip.label)),
#'                                  vector.Theta = vector.Model_Simulation_Theta)
#'
#'
#' handle.Sim_Results <- Function_ConductSimExperimentOU_ModelTestBM_OU(list.Model_Simulation_OU = list.Model_Simulation_OU, numeric.NumberOfReps = 5)
#'
#' ###################
#' # Compare this to #
#' ###################
#' ##########################
#' # Build simulation model #
#' ##########################
#' phylogeny.RandomTree <- pbtree(b = 1, d = 0, n = 100)
#' vector.Model_Simulation_Theta <- c(1, 0.0001)
#' names(vector.Model_Simulation_Theta) <- c("Sig2", "alpha")
#'
#' list.Model_Simulation_OU <- list(handle.Phylogeny = phylogeny.RandomTree,
#'                                  string.Model = "OU",
#'                                  vector.Z = rep(0, length(phylogeny.RandomTree$tip.label)),
#'                                  vector.Theta = vector.Model_Simulation_Theta)
#'
#'
#' handle.Sim_Results <- Function_ConductSimExperimentOU_ModelTestBM_OU(list.Model_Simulation_OU = list.Model_Simulation_OU, numeric.NumberOfReps = 5)


##################################################
# Function_ConductSimExperimentOU_ModelTestBM_OU #
##################################################
Function_ConductSimExperimentOU_ModelTestBM_OU <- function(list.Model_Simulation_OU, numeric.NumberOfReps){

  ###################
  # Summarize input #
  ###################
  matrix.SimulationResults <- matrix(nrow = numeric.NumberOfReps, ncol = 4)
  colnames(matrix.SimulationResults) <- c("P-value", "LnLRatio", "dH", "dKL")

  ##########################
  # Build Simulation Model #
  ##########################
  handle.Model_Simulation_Tree <- list.Model_Simulation_OU$handle.Phylogeny
  string.Model_Simulation_Model <- list.Model_Simulation_OU$string.Model
  vector.Model_Simulation_Z <- list.Model_Simulation_OU$vector.Z
  vector.Model_Simulation_Theta <- list.Model_Simulation_OU$vector.Theta
  numeric.Model_Simulation_Sig2 <- vector.Model_Simulation_Theta["Sig2"]
  numeric.Model_Simulation_alpha <- vector.Model_Simulation_Theta["alpha"]

  ##############################
  # Transform phylogeny to VCV #
  ##############################
  phylogeny.Rescaled_Simulation_Model <- rescale(x = handle.Model_Simulation_Tree, model = "OU")
  phylogeny.Rescaled_Simulation_Model <- phylogeny.Rescaled_Simulation_Model(alpha = numeric.Model_Simulation_alpha, sigsq = numeric.Model_Simulation_Sig2)
  matrix.VCV_Rescale_Simulation_Model <- vcv(phy = phylogeny.Rescaled_Simulation_Model)

  #######################
  # Conduct simulations #
  #######################
  for (i in 1:numeric.NumberOfReps){

    #################
    # Simulate data #
    #################
    #vector.Data_Simulated_CharacterTraits <- fastBM(tree = phylogeny.Rescaled_Simulation_Model, a = 0, mu = 0)
    vector.Data_Simulated_CharacterTraits <- mvrnorm(n = 1, mu = vector.Model_Simulation_Z, Sigma = matrix.VCV_Rescale_Simulation_Model, tol = 1e-6)

    ##########################################
    # Estimate a BM model for simulated data #
    ##########################################
    handle.ModelFit_BM <- fitContinuous(phy = handle.Model_Simulation_Tree, dat = vector.Data_Simulated_CharacterTraits, SE = 0, model = "BM")
    handle.ModelFit_OU <- fitContinuous(phy = handle.Model_Simulation_Tree, dat = vector.Data_Simulated_CharacterTraits, SE = 0, model = "OU", bounds = list(alpha=c(10^-10,100)))

    ###############################
    # Build estimate for BM model #
    ###############################
    numeric.Model_BM_Sig2 <- handle.ModelFit_BM$opt$sigsq
    numeric.Model_BM_z0 <- handle.ModelFit_BM$opt$z0
    vector.Model_BM_Theta <- c(numeric.Model_BM_Sig2)
    names(vector.Model_BM_Theta) <- c("Sig2")
    list.Model_BM <- list(handle.Phylogeny = handle.Model_Simulation_Tree,
                          string.Model = "BM",
                          vector.Z = rep(numeric.Model_BM_z0, length(handle.Model_Simulation_Tree$tip.label)),
                          vector.Theta = vector.Model_BM_Theta)


    ###############################
    # Build estimate for OU model #
    ###############################
    numeric.Model_OU_Sig2 <- handle.ModelFit_OU$opt$sigsq
    numeric.Model_OU_alpha <- handle.ModelFit_OU$opt$alpha
    numeric.Model_OU_z0 <- handle.ModelFit_OU$opt$z0
    vector.Model_OU_Theta <- c(numeric.Model_OU_Sig2, numeric.Model_OU_alpha)
    names(vector.Model_OU_Theta) <- c("Sig2", "alpha")
    list.Model_OU <- list(handle.Phylogeny = handle.Model_Simulation_Tree,
                          string.Model = "OU",
                          vector.Z = rep(numeric.Model_OU_z0, length(handle.Model_Simulation_Tree$tip.label)),
                          vector.Theta = vector.Model_OU_Theta)

    ####################
    # Compute distance #
    ####################
    vector.Results_Distance_Between_EstimatedModels <- Function_ComputeDistances(list.Model_01 = list.Model_BM, list.Model_02 = list.Model_OU)
    matrix.SimulationResults[i, "dH"] <- vector.Results_Distance_Between_EstimatedModels["dH"]
    matrix.SimulationResults[i, "dKL"] <- vector.Results_Distance_Between_EstimatedModels["dKL"]

    #####################
    # Compute LnL ratio #
    #####################
    numeric.LnL_Given_BM <- handle.ModelFit_BM$opt$lnL
    numeric.LnL_Given_OU <- handle.ModelFit_OU$opt$lnL
    numeric.LnL_Ratio <- 2*(numeric.LnL_Given_OU-numeric.LnL_Given_BM)
    matrix.SimulationResults[i, "LnLRatio"] <- numeric.LnL_Ratio

    #################
    # Compute chisq #
    #################
    numeric.Pvalue_ChiSq <- pchisq(numeric.LnL_Ratio,1,lower.tail=F)
    matrix.SimulationResults[i,"P-value"] <- numeric.Pvalue_ChiSq

  }
  return(matrix.SimulationResults)
}
