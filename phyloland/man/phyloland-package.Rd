\name{phyloland-package}
\alias{phyloland-package}
\alias{phyloland}
\docType{package}

\title{
  Modelling competition and  dispersal in a statistical phylogeographic framework
}

\description{
The R package "phyloland" implements a landscape colonization model in a statistical phylogeographic framework, explicitely modelling dispersal and competition. 
The model aims at testing whether these two factors have had a significant impact on the process of landscape colonization that explains the current geographical distribution of the evolutionary units under study.
Given (i) their phylogeny and (ii) their geographic location, the distributions of these two parameters is estimated.

The landscape is defined in a 2-dimensional space (latitude and longitude) as a set of locations that can be colonized through a sequence of dated dispersal events.
According to this approach, dispersals are mapped on the genealogy, each internal node corresponding to a specific dispersal event.

This model takes into account two constraints on dispersals : 
\itemize{
  \item{The geographical distance separating the source and the post-dispersal locations:}{ does the dispersal probability decrease with the geographic distance separating the locations of the landscape?
  A normal distribution centered on the source location is used to model dispersal.
  The standard deviation sigma determines the geographical dispersal bias.
  Low values for sigma indicate that the units of interest tend to disperse locally.
  Conversely, with highly mobile species or limited space, the probability to migrate to any locations in the landscape may be uniform, corresponding to larger values of sigma.}
  \item{The occupancy state of the post-dispersal location:}{ does the probability of colonizing a location vary with its occupancy state, i.e. whether it has previously been colonized or not?
  An occupied location is colonized with probability lambda.
  If competition for space takes place, dispersals to occupied locations are prevented and the estimates of the parameter lambda will fall below 1.}
}

The 3 model parameters are : 
\itemize{
 \item{sigma}{: the standard deviation of the normal distribution for dispersal kernel.}
 \item{lambda}{: the probability of migrating to an occupied location.}
 \item{tau}{: the overall dispersal rate, the number of dispersals per unit of time.}
}

Bayesian Markov chain Monte Carlo (MCMC) algorithms are used to sample the model parameters and the geographic locations at the internal ndoes of the phylogeny.
Uniform distributions are used as prior distribution for all parameters.
}

\details{
The following package functions allow the user to implement the model and investigate the results : 
\itemize{
  \item{The function \link{PLD_interface} performs the Bayesian estimations from a phylogeny and a location file.}

  \item{\link{PLD_plot_trees} plots a number of sampled trees with their internal and tips locations.}

  \item{\link{PLD_loc_mrca} returns the posterior distribution of each location of the Most Recent Common Ancestor of a set of tips. A barplot can be displayed.}

  \item{\link{PLD_stat_mig} lists the dispersal events from the genealogies sampled by the MCMC.}

  \item{\link{PLD_plot_stat_mig} displays the density plots of the dispersal time points for each location of the space.}
}
}

