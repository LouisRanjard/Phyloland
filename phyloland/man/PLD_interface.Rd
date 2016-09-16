\name{PLD_interface}

\alias{PLD_interface}

\title{
Model parameters and internal locations estimations
}

\description{
Function that estimates model parameters, genealogies and internal locations through Bayesian Markov chain Monte Carlo (MCMC) algorithm.
}

\usage{
PLD_interface(fileTREES, fileDATA, num_step = 1e+05, freq = 100, 
burnin = 0, ess_lim = 100, sigma = NA, lambda = NA, tau = NA, 
num_step_sigma = 1, num_step_lambda = 1, num_step_tau = 1, 
id_filena = NA, pattern_trees_likelihood = "treeLikelihood", 
names_locations = NA)
}

\arguments{
  \item{fileTREES}{
a character string that specifies the name of the file containing the phylogenetic trees.
A Nexus file (e.g. BEAST output) with tips names (no space characters allowed), the phylogenetic trees and the trees likelihood. 
If only one tree is contained in the file, the genealogies are not sampled in the MCMC.
}

  \item{fileDATA}{
a character string that specifies the name of the file containing the tips locations : a text file with 3 columns (separated by tabs).
The first one with the tips names (the same as in \code{fileTrees}), the second one with the location latitudes and the last one with the location longitudes (both in decimal degrees).
No header.
}

  \item{num_step}{
a strictly positive integer specifying the length of the Markov Chain (the number of MCMC steps).
Suggested values : 50 000 - 100 000.
}

  \item{freq}{
a strictly positive integer specifying how often the sampled model parameters and the locations are saved.
Suggested values : 50 - 100.
}

  \item{burnin}{
a strictly positive integer specifying the number of trees dropped from the file \code{fileTREES}.
Suggested value : 10\% of the total number of trees, 0 in case of only one tree.
}

  \item{ess_lim}{
a strictly positive integer specifying the Effective Sample Size.
The MCMC stops when the ess of all parameters reach this value.
}

  \item{sigma}{
a vector of two elements, strictly positive values specifying the standard deviation of the normal distribution for dispersal steps in each dimension.
Default is NA: the 2 sigma parameters are estimated in the MCMC.
Alternatively, if values are specified the sigma parameters are not estimated by the MCMC.
}

  \item{lambda}{
a strictly positive value specifying the probability of migrating to an occupied location.
Default is NA: the lambda parameter is estimated in the MCMC.
Alternatively, if a value is specified the lambda parameter is not estimated by the MCMC.
A lambda value equal to 1 indicates that competitive exclusion did not impact the colonization process.
}

  \item{tau}{
a strictly positive value specifying the overall dispersal rate.
Default is NA: the tau parameter is estimated in the MCMC.
Alternatively, if a value specified the tau parameter is not estimated by the MCMC.
}

  \item{num_step_sigma}{
a strictly positive integer specifying how many times the sigma parameters are sampled at each MCMC step.
}

  \item{num_step_lambda}{
a strictly positive integer specifying how many times the lambda parameter is sampled at each MCMC step.
}

  \item{num_step_tau}{
a strictly positive integer specifying how many times the tau parameter is sampled at each MCMC step.
}

  \item{id_filena}{
an optional character string specifying the ID for the output files.
}

  \item{pattern_trees_likelihood }{
a character string specifying the text pattern to retrieve the trees likelihood in \code{fileTREES}.
}

  \item{names_locations}{
a vector of character strings that specifies the location names, in the same order as in \code{fileDATA}.  
}
}

\details{
MCMC runs stop when the \code{num_step} is reached or when the ess of all estimated parameters equal \code{ess_lim}.

The respective ess for estimating convergence of the chains regarding each parameter is calculated by monitoring the autocorrelation between successive sampled states in the MCMC.

%Six files are written : the first one with the trees sampled (txt file), the second one with the locations sampled (log file), the third one with the model parameters sampled (log file), the fourth one with the simulation parameters (csv file), the fifth one with informations about the model parameters posterior distributions (csv file) and the sixth one with the object "phyloland" (Rdata file).
}

\value{
an object of class "phyloland" with the following components:

\item{trees}{  an object of class "multiPhylo" or "phylo" (package ape) containing the trees sampled in the MCMC.
If the \code{ess_lim} is not reached, the number of trees is \code{num_step}/\code{freq}.}

\item{locations}{  a matrix containing the internal locations sampled in the MCMC and the tips locations (columns) for each sampled tree(rows).}

\item{tips}{  a vector of character strings containing the tips names (from \code{fileTREES} and \code{fileDATA}).}

\item{space}{  a matrix containing the unique locations from \code{fileDATA}.}

\item{mcmc}{  a list containing the posterior distributions of each parameters.}

\item{sigma_limit}{  a vector indicating the threshold for sigma below which limited dispersal occurs.}
}

