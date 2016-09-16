setwd("~/Documents/phylogeo/SVNphyloland/trunk/phyloland_tutorial/")
library(phyloland)

# Location names for the tips in the tree (same order as in tree file, tree_Banza.nex)
names_locations = c("Maui_Nui", "Maui_Nui", "Maui_Nui", "Maui_Nui", "Kauai", "Kauai", "Maui_Nui", "Maui_Nui", "Maui_Nui", "Maui_Nui", "Nihoa", "Nihoa", "Hawaii", "Hawaii", "Hawaii", "Oahu", "Oahu", "Maui_Nui", "Maui_Nui", "Oahu", "Oahu")

# Run the model, on one consensus tree:
#Banza = PLD_interface(fileTREES="tree_Banza.nex", fileDATA="locations_Banza.txt", num_step=1e5, freq=1e2, ess_lim=5e2, names_locations=names_locations)
# or alternatively on the Beast posterior distribution of trees:
Banzap = PLD_interface(fileTREES="tree_Banza_posterior.nex", fileDATA="locations_Banza.txt", num_step=1e3, freq=1e2, ess_lim=5e2, names_locations=names_locations)


# Plot sampled internal node locations
#pdf("figures/tree_20_1000_Banza.pdf", height=6, width=12)
#PLD_plot_trees(x = Banza, sub_sample = 1)
PLD_plot_trees(x = Banza, sub_sample = c(20, 1000), one_plot = TRUE)
#dev.off();

# Plot estimated parameter values
#mcmc_log = read.csv("tracer_Banza.log",sep="\t")
#mcmc_summary = read.csv("phyloland_Banza.csv",sep=",",row.names=NULL)
mcmc_log = Banza$mcmc[[1]]
x11(); 
#pdf("figures/post_dist.pdf", height=6, width=6);
par(mfrow=c(2,2))
hist(mcmc_log[,2],40,xlab="sigma 1",main="sigma 1"); abline(v=Banza$mcmc[[8]][1],lwd=2,col="red")
hist(mcmc_log[,3],40,xlab="sigma 2",main="sigma 2"); abline(v=Banza$mcmc[[8]][2],lwd=2,col="red")
hist(log(mcmc_log[,4]),40,xlab="log(lambda)",main="lambda"); abline(v=0,lwd=2,col="red")
hist(mcmc_log[,5],40,xlab="tau",main="tau");
#dev.off();

# Plot Most Recent Common Ancestor locations
locations = PLD_loc_mrca(x = Banza, tips = Banza$tips, sub_sample = 0, plot_distrib = TRUE)
locations$frequencies

# Plot the oldest migration time to each location
stat = PLD_stat_mig(x = Banza, sub_sample = 1:100, first = TRUE)
x11(); 
#pdf("figures/timemig_Banzap.pdf", height=6, width=6);
PLD_plot_stat_mig(timemat = stat$timemat, color = NA, xy_legend = NA)
#dev.off();

