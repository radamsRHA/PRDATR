% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/Function_ComputePairwiseDistances_SimulateTips_Depth.R
\name{Function_ComputePairwiseDistances_SimulateTips_Depth}
\alias{Function_ComputePairwiseDistances_SimulateTips_Depth}
\title{Function_ComputePairwiseDistances_SimulateTips_Depth: function to generate pairwise distance computations}
\usage{
Function_ComputePairwiseDistances_SimulateTips_Depth(
  numeric.NumberReps,
  numeric.n,
  numeric.birth,
  string.Model,
  vector.Theta
)
}
\arguments{
\item{numeric.NumberReps}{Number of simulation replicates}

\item{numeric.n}{Number of tips}

\item{numeric.birth}{Birth Rate}

\item{string.Model}{String for evolutonary model}
}
\value{
List list with (1) matrix.PairwiseDistances_H, (2) matrix.PairwiseDistances_KL, and (3)vector.PairwiseDistances_H
}
\description{
This function returns a list with (1) matrix.PairwiseDistances_H, (2) matrix.PairwiseDistances_KL, and (3)vector.PairwiseDistances_H
}
\examples{



################################
# Simulate datasets with n = 3 #
################################
vector.Theta <- c(1, 1)
names(vector.Theta) <- c("Sig2", "alpha")
list.RESULTS_n3_b1_a1 <- Function_ComputePairwiseDistances_SimulateTips(numeric.NumberReps = 100,
                                                                        numeric.n = 3, numeric.birth = 1,
                                                                        string.Model = "OU",
                                                                        vector.Theta = vector.Theta)$vector.PairwiseDistances_H
}
\keyword{brownian}
\keyword{continuous}
\keyword{distance,}
\keyword{distances,}
\keyword{evolution}
\keyword{model}
\keyword{motion,}
\keyword{phylogenetic}
\keyword{probabilistic}
\keyword{trait}
