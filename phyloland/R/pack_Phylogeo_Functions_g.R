#Function that estimates model parameters, genealogies and internal locations through Bayesian Markov chain Monte Carlo (MCMC) algorithms.
#fileTREES: contain tree in nexus format, tip names must be the same as in fileDATA
#fileDATA: for each tip name, give the corresponding space coordinates
PLD_interface <- function(fileTREES, fileDATA, num_step=100000, freq=100, burnin=0, ess_lim=100, sigma=NA, lambda=NA, tau=NA, num_step_sigma=1, num_step_lambda=1, num_step_tau=1, id_filena=NA, pattern_trees_likelihood="treeLikelihood", names_locations=NA){
  if(missing(fileTREES)){
    stop("'fileTREES' is missing")	
  }
  if(missing(fileDATA)){
    stop("'fileDATA' is missing")	
  }
  verif_arg(fileTREES, type = "character", len = 1, name_function = "PLD_interface")
  verif_arg(fileDATA, type = "character", len = 1, name_function = "PLD_interface")
  if (sum(!(is.na(sigma))) != 0){verif_arg(sigma, type = "numeric", len = 2, name_function = "PLD_interface", positive = 2)}
  if (sum(!(is.na(lambda))) != 0){verif_arg(lambda, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)}
  if (sum(!(is.na(tau))) != 0){verif_arg(tau, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)}
  verif_arg(burnin, type = "numeric", len = 1, name_function = "PLD_interface", positive = 1)
  verif_arg(num_step, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(freq, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(num_step_sigma, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(num_step_lambda, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(num_step_tau, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  verif_arg(ess_lim, type = "numeric", len = 1, name_function = "PLD_interface", positive = 2)
  if (sum(!(is.na(id_filena))) != 0){verif_arg(id_filena, type = "character", len = 1, name_function = "PLD_interface")}
  verif_arg(pattern_trees_likelihood, type = "character", len = 1, name_function = "PLD_interface")
  
  for(i in 1:length(readLines(fileDATA))){
    if (length(strsplit(readLines(fileDATA)[i],split="\t")[[1]]) != 3){
      stop("Error in 'fileDATA', is it TAB delimited?")	
    }
  } 
  if (sum(is.na(as.matrix(read.table(fileDATA, header = FALSE, sep="\t"))))){
    stop("Error in 'fileDATA', is it TAB delimited?")	
  }
  
  # read trees phylo
  trees_phylo = read.nexus(fileTREES)
  if (is(trees_phylo,"phylo")){
    sample_geneal = c(trees_phylo)
    if (sum(trees_phylo$edge.length<0)>0) {
      stop('[Phyloland:PLD_interface] non-clock-like phylogenetic tree (negative branch length)')
    }
    #gtreelikelihood = 1
  }else if (is(trees_phylo,"multiPhylo")){
    sample_geneal = trees_phylo[(burnin+1):length(trees_phylo)]
    #trees_par = readLines(fileTREES)[which(regexpr("STATE",readLines(fileTREES)) != -1)]
    #pattern = paste(pattern_trees_likelihood,"=",sep = "")
    # tree likelihood (in BEAST file)
    #gtreelikelihood = lapply(trees_par, treeLikeli, pattern = pattern)
  }
  
  #tips, watch out: because of the way ape read nexus file, the order of tips may be different in each tree
  #unless there is a "Translate" block in nexus file
  tips = sample_geneal[[1]]$tip.label
  space_loc = read_space(fileDATA, tips)	#read space file
  loc_tips = space_loc[[1]] #give for each element of tips the index of location in space
  space = space_loc[[3]]
  
  if (sum(is.na(names_locations)) != 0){
    #names_locations = as.character(c(1:dim(space)[2]))
    names_locations = as.character(tips)
  }else{
    names_locations = as.character(names_locations)
  }
  verif_arg(names_locations, len = length(tips), name_function = "PLD_interface")
  colnames_space = rep("",dim(space)[2])
  for (n in 1:dim(space)[2]){
    colnames_space[n] = names_locations[which(loc_tips==n)[1]]
  }
  colnames(space) = colnames_space
  
  # gtreel and locations_sim for mcmc_phyloland
  gtreel = converttree(sample(sample_geneal, 1)[[1]])
  locations_sim = cbind(rep(0,dim(gtreel$nodes)[1]),rep(0,dim(gtreel$nodes)[1]),rep(0,dim(gtreel$nodes)[1]))
  locations_sim[1:length(loc_tips),1] = loc_tips
  
  # files creation
  if (is.na(id_filena)){
    #id_filena = format(Sys.time(),"%a%d%b%Y_%H%M%S")
    filena = NA
    filena_param = NA
    filena_loc = NA
    filena_tracer = NA
    filena_tree = NA
  }else{
    filena = paste('phyloland_',id_filena,'.csv',sep='')
    filena_param = paste('parameters_',id_filena,'.csv',sep='')
    filena_loc = paste('loc_',id_filena,'.log',sep='')
    filena_tracer = paste("tracer_",id_filena,".log",sep="")
    filena_tree = paste("tree_",id_filena,".tre",sep="")
  }
  
  # no output filename given, should not be allowed...
  if (!is.na(filena_loc)){
    nodes = c(1:length(gtreel$nodes[,1]))
    nodes[1:length(tips)] = tips
    write.table(t(c("num_step",nodes)),row.names=FALSE,col.names=FALSE,file=filena_loc,sep="\t",eol="\n",append=FALSE)
  }
  
  #load tree likelihood
  treelikelihood = read_treelikelihood(filena_tree)

  #dmethod="euclidean"
  dmethod="distkm"  # get dispersal kernel in km
  
  # MCMC
  mcmco = mcmc_phyloland(space=space, gtreel=gtreel, simul_values=list(sigma,lambda,tau,locations_sim), treelikelihood, est_sigma=c(1,1), est_lambda=1, est_tau=1, sample_geneal=sample_geneal, sample_loc=1, plot_phylo=0, plot_MCMC=0, save_MCMC=0, Nstep=num_step, freq=freq, Nstep_sigma=num_step_sigma, Nstep_lambda=num_step_lambda, Nstep_loc=1, Nstep_genealogy=1, Nstep_tau=num_step_tau, pchange=0, ess_lim=ess_lim, model=4, show_loc=0, filena_loc, filena_tracer, filena_tree, dmethod)
  
  stat_res = stat_phyloland(mcmco[[1]])
  stat_loc = vector("numeric",7)
  
  trees = mcmco[[3]]
  #sampled_loc = read_loc(filena_loc)
  sampled_loc = mcmco[[2]][,2:ncol(mcmco[[2]])]
  ind = mcmco[[6]]
  rownames(space) = c("Dimension1", "Dimension2")
  x <- (list( trees = trees, locations = sampled_loc, tips = tips, space = space, index = ind, mcmc = mcmco, sigma_limit = round(unlist(mcmco[[8]]),3)))
  class(x) <- "phyloland"	
  # phyloland object
  #save(x, file = paste("phyloland_",id_filena,".Rdata",sep = ""))
  
  # print dispersal kernels
  if (dmethod=="distkm") {
    #print(stat_res)
    #browser()
    #if ( as.numeric(unlist(strsplit(stat_res,","))[9])<mcmco[[8]][1] ) {
      print(paste("Latitudinal dispersal kernel (50% quantile):",format(round(qnorm(0.75,mean=0,sd=sqrt(as.numeric(unlist(strsplit(stat_res,","))[9]))), 2), nsmall=2),"km"))
    #}
    #if ( as.numeric(unlist(strsplit(stat_res,","))[19])<mcmco[[8]][2] ) {
      print(paste("Longitudinal dispersal kernel (50% quantile):",format(round(qnorm(0.75,mean=0,sd=sqrt(as.numeric(unlist(strsplit(stat_res,","))[19]))), 2), nsmall=2),"km"))
    #}
  }
  
  # files
  if (!is.na(filena_param)){
    cat(paste('num_location',sep=","),file=filena_param,sep=",",append=FALSE)
    for (ndim in 1:dim(space)[1]){
      cat(paste(',prior_sigma',ndim,sep=""),file=filena_param,sep=",",append=TRUE)
    }
    cat(paste(',prior_lambda', 'prior_tau', 'num_step', 'freq', 'num_step_sigma', 'num_step_lambda', 'num_step_tau', 'ess_lim', '\n', sep=","),file=filena_param,sep=",",append=TRUE)
    
    cat(dim(space)[2],file=filena_param,sep=",",append=TRUE)
    
    if (sum(is.na(sigma)) != 0){
      for (ndim in 1:dim(space)[1]){
        cat(paste(',unif',sep=","),file=filena_param,sep=",",append=TRUE)
      }
    }else{
      cat('',sigma,file=filena_param,sep=",",append=TRUE)
    }
    if (is.na(lambda)){
      cat(paste(',unif',sep=","),file=filena_param,sep=",",append=TRUE)
    }else{
      cat('',lambda,file=filena_param,sep=",",append=TRUE)
    }	
    if (is.na(tau)){
      cat(paste(',unif',sep=","),file=filena_param,sep=",",append=TRUE)
    }else{
      cat('',tau,file=filena_param,sep=",",append=TRUE)
    }		
    cat(c('',num_step, freq, num_step_sigma, num_step_lambda, num_step_tau, ess_lim,'\n'),file=filena_param,sep=",",append=TRUE)
  }
  
  if (!is.na(id_filena)){
    for (ndim in 1:dim(space)[1]){
      if (ndim == 1){
        cat(paste(paste('min_sigma',ndim,sep=""),'quant05','quant25','quant50','quant75','quant95','max','mode','mean','last',sep=","),file=filena,sep=",",append=FALSE)
      }else{
        cat(paste(paste(',min_sigma',ndim,sep=""),'quant05','quant25','quant50','quant75','quant95','max','mode','mean','last',sep=","),file=filena,sep=",",append=TRUE)
      }
    }
    cat(paste(',min_lambda','quant05','quant25','quant50','quant75','quant95','max_lambda','mode_lambda','mean_lambda','last_lambda',sep=","),file=filena,sep=",",append=TRUE)
    cat(paste(',min_tau','quant05','quant25','quant50','quant75','quant95','max_tau','mode_tau','mean_tau','last_tau',sep=","),file=filena,sep=",",append=TRUE)
    for (ndim in 1:dim(space)[1]){
      cat(paste(',sigma_upper_limit',ndim,sep=""), file=filena, sep=",", append=TRUE)
    }
    cat('\n',file=filena, sep=",",append=TRUE)
    cat(c(unlist(stat_res),round(unlist(mcmco[[8]]),3),'\n'),file=filena,sep=",",append=TRUE)
  }
  
  return(x)
}


# Function that plots sampled trees with their locations.
PLD_plot_trees <- function(x, sub_sample = 0, one_plot = FALSE){
  if (missing(x)){
    stop("'x' is missing")	
  }
  verif_arg(x, type = "phyloland", name_function = "plot_trees")
  if(class(x$trees) == "phylo"){
    x$trees = c(x$trees)	
  }	
  verif_arg(sub_sample, type = c("numeric","integer"), name_function = "plot_trees")
  verif_arg(one_plot, type = "logical", len = 1, name_function = "plot_trees")
  trees = x$trees
  verif_arg(x$trees, type = "multiPhylo" , name_function = "plot_trees")
  sampled_loc = x$locations
  
  if(length(sub_sample) == 1){
    if (sub_sample == 0){
      sub_sample = 1:length(trees)	
    }		
  }	
  if (sum(is.element(sub_sample,1:length(trees))==FALSE) != 0){
    stop(paste(sub_sample[which(is.element(sub_sample,1:length(trees))==FALSE)[1]]," : subscript out of bounds"))
  }
  
  s = sub_sample
  
  n = length(s)
  if (one_plot == TRUE ){
    nb_row = round(sqrt(n))
    nb_col = ceiling(sqrt(n))
    par(mfrow = c(nb_row,nb_col))
  }
  for (i in 1:n){
    print(paste("tree",s[i]),sep = "")
    plot_trees_tips(trees[[s[i]]],sampled_loc[s[i],],x$space)
    if(one_plot == FALSE){
      if (i < n){
        readline()
      }
    }
  } 
}

# Function that displays the density plots of the migration times for each location.
PLD_loc_mrca <- function(x, tips, sub_sample = 0, plot_distrib = FALSE, col = NA){
  if (missing(x)){
    stop("'x' is missing")	
  }
  verif_arg(x, type = c("phyloland"), name_function = "loc_mrca")
  if(class(x$trees) == "phylo"){
    x$trees = c(x$trees)	
  }
  trees = x$trees
  verif_arg(x$trees, type = "multiPhylo" , name_function = "loc_mrca")
  sampled_loc = x$locations
  corresp = x$tips
  space = x$space
  names_locations = colnames(space)
  verif_arg(plot_distrib, type = "logical", len = 1, name_function = "loc_mrca")
  verif_arg(sub_sample, type = c("numeric","integer"), name_function = "loc_mrca")
  if(length(sub_sample) == 1){
    if (sub_sample == 0){
      sub_sample = 1:length(trees)	
    }		
  }
  if (sum(is.element(sub_sample,1:length(trees))==FALSE) != 0){
    stop(paste(sub_sample[which(is.element(sub_sample,1:length(trees))==FALSE)[1]]," : subscript out of bounds"))
  }
  trees = trees[sub_sample]
  sampled_loc = sampled_loc[sub_sample,]
  
  if (missing(tips)){
    tips = corresp
  }else{
    tips = unique(tips)	
  }
  verif_arg(tips, len = c(1:length(corresp)), name_function = "loc_mrca")
  tips_ind = vector("numeric",length(tips))
  tips = as.character(tips)	
  for (i in 1:length(tips)){
    if (length(which(corresp == tips[i]))==0){
      stop(paste(" tip \"" , tips[i],"\" not found", sep = ""))
    }else{
      tips_ind[i] = which(corresp == tips[i])
    }
  }	
  if (sum(is.na(col))!=0){
    col = NULL	
  }
  loc_MRCA = vector("numeric",length(trees))	
  if (length(tips) == 1){
    mat = c(space[1,sampled_loc[1,tips_ind]], space[2,sampled_loc[1,tips_ind]], 1)
    loc_MRCA = rep(sampled_loc[1,tips_ind],length(trees))
    res = list(frequencies = mat , locationsMRCA = loc_MRCA)
    if (plot_distrib == TRUE){
      tab <- table(loc_MRCA)/sum(table(loc_MRCA))
      tab_complete <- vector("numeric",dim(space)[2])
      tab_complete[as.integer(rownames(tab))] = tab
      names(tab_complete) = names_locations
      #x11() #dev.new()
      barplot(tab_complete, main = "", ylab = "frequencies", xlab = "locations" , col = col)
    }
    return(res)	
  }else{
    for (i in 1:length(trees)){
      loc_MRCA[i] = loc_ancestor(converttree(trees[[i]]), tips_ind, sampled_loc[i,])
    }
    mat <- matrix(ncol = 3, nrow = length(table(loc_MRCA)), dimnames = list(as.character(names(sort(table(loc_MRCA),decreasing=TRUE))),c(rownames(space),"Frequency")))
    sort_loc = as.integer(names(sort(table(loc_MRCA),decreasing=TRUE)))
    mat[,1] = space[1,sort_loc]
    mat[,2] = space[2,sort_loc]
    mat[,3] = sort(table(loc_MRCA),decreasing=TRUE) / sum(table(loc_MRCA))
    rownames(mat) = names_locations[as.integer(rownames(mat))]
    if (plot_distrib == TRUE){				
      tab <- table(loc_MRCA)/sum(table(loc_MRCA))
      tab_complete <- vector("numeric",dim(space)[2])
      tab_complete[as.integer(rownames(tab))] = tab
      names(tab_complete) = names_locations
      #x11() #dev.new()
      barplot(tab_complete, ylab = "frequencies", xlab = "locations" , col = col)
    }
    #if (plot_distrib[2] == TRUE){
    #			x11()
    #			loc_lat = vector("numeric", length(loc_MRCA))
    #			loc_long = vector("numeric", length(loc_MRCA))
    #			for(i in 1:length(loc_MRCA)){
    #				loc_lat[i] = space[1,loc_MRCA[i]]
    #				loc_long[i] = space[2,loc_MRCA[i]]
    #			}
    #			kde = kde2d(loc_lat, loc_long, n = 300)
    #			image(kde)
    #			points(space[1,],space[2,],cex=2, pch = 16)
    #		}
    res = list(frequencies = mat, locationsMRCA = loc_MRCA)
    return(res)
  }
}

# Function that lists the migration events from the genealogies sampled by the MCMC.
PLD_stat_mig <- function(x, sub_sample = 0, first = FALSE){
  if (missing(x)){
    stop("'x' is missing")	
  }
  verif_arg(x, type = "phyloland", name_function = "stat_mig")
  if(class(x$trees) == "phylo"){
    x$trees = c(x$trees)	
  }
  verif_arg(sub_sample, type = c("numeric","integer"), name_function = "stat_mig")
  verif_arg(first, type = "logical", len = 1, name_function = "stat_mig")
  trees = x$trees
  verif_arg(x$trees, type = "multiPhylo" , name_function = "stat_mig")
  sampled_loc = x$locations
  space = x$space
  names_locations = colnames(space)
  
  if(length(sub_sample) == 1){
    if (sub_sample == 0){
      sub_sample = 1:length(trees)	
    }		
  }
  if (sum(is.element(sub_sample,1:length(trees))==FALSE) != 0){
    stop(paste(sub_sample[which(is.element(sub_sample,1:length(trees))==FALSE)[1]]," : subscript out of bounds"))
  }
  trees = trees[sub_sample]
  sampled_loc = sampled_loc[sub_sample,]
  
  migmat = matrix(0,nrow=dim(space)[2],ncol=dim(space)[2])
  if (first == FALSE){
    timemat = array(0,c(dim(space)[2],dim(space)[2],(sum(is.na(converttree(trees[[1]])[,2]))-1)*length(trees)))
  }else{
    timemat = array(0,c(dim(space)[2],dim(space)[2],length(trees)))
  }
  for (m in 1:length(sampled_loc[,1])){
    treel = reorder_treel(converttree(trees[[m]]))[[1]]
    history = treel2mig(treel,sampled_loc[m,],space)
    # record location at the root = ancestral location of MRCA
    timemat[history[1,1],history[1,1],m] = history[1,5+dim(space)[1]] # ! will be confounded with the first migration event
    for (n in 1:dim(history)[1]){
      migmat[history[n,1],history[n,2]] = migmat[history[n,1],history[n,2]] + 1
      if (first == FALSE){
        timemat[history[n,1],history[n,2],which(timemat[history[n,1],history[n,2],]==0)[1]] = history[n,5+dim(space)[1]]
      }else{
        if (sum(timemat[,history[n,2],m])==0){
          timemat[history[n,1],history[n,2],m] = history[n,5+dim(space)[1]]
        }
      }
    }
  }
  rownames(timemat) = names_locations
  colnames(timemat) = names_locations
  rownames(migmat) = names_locations
  colnames(migmat) = names_locations
  statmig = list(migmat = migmat, timemat = timemat)
  return(statmig)
}

#Function that displays the density plots of the migration times for each location.
PLD_plot_stat_mig <- function(timemat, color = NA, xy_legend = NA, group = 0){
  if (missing(timemat)){
    stop("'timemat' is missing")	
  }
  verif_arg(timemat, type = c("array"), name_function = "plot_stat_mig")
  if (sum(is.na(xy_legend)) == 0){
    verif_arg(xy_legend, type = "numeric", len = 2, name_function = "plot_stat_mig")
  }
  if (sum(is.na(color)) == 0){
    verif_arg(color, len = nrow(timemat), name_function = "plot_stat_mig")
  }
  if (group==1) {
    # Group locations according to row names
    timemat2 = array(0,c(length(unique(rownames(timemat))),length(unique(rownames(timemat))),length(timemat[1,1,])))
    rownames(timemat2) = unique(rownames(timemat))
    colnames(timemat2) = unique(colnames(timemat))
    for (n in 1:length(timemat[1,1,])){
      for (loc1 in unique(rownames(timemat))){
        for (loc2 in unique(colnames(timemat))){
          timemat2[loc1,loc2,n] = max( timemat[which(rownames(timemat[,,1])==loc1),which(colnames(timemat[,,1])==loc2),n] )
        }
      }
    }
    timemat = timemat2
  }
  junk.x = NULL
  junk.y = NULL
  for(i in 1:(nrow(timemat))){
    if (length(timemat[,i,][timemat[,i,]!=0])>0){
      junk.x = c(junk.x, density(timemat[,i,][timemat[,i,]!=0])$x)
      junk.y = c(junk.y, density(timemat[,i,][timemat[,i,]!=0])$y)
    }
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  if (sum(is.na(color)) != 0){
    color = rainbow(nrow(timemat))
  }	
  tmp = density(timemat[,1,][timemat[,1,]!=0])
  ymax = max(tmp$y)
  plot(tmp, col = color[1], xlim = xr, ylim = yr, main = "", ylab = "Density", xlab ="Time (BP)", lwd = 2)
  if (nrow(timemat) > 1){
    for(i in 2:nrow(timemat)){
      if (length(timemat[,i,][timemat[,i,]!=0])>0){
        tmp = density(timemat[,i,][timemat[,i,]!=0])
        if (max(tmp$y)>ymax) ymax = max(tmp$y)
        lines(tmp, col = color[i], lwd = 2)
      }
    }
  }
  if (sum(is.na(xy_legend))!=0){
    xy_legend = c(0.75*max(xr),max(yr))	
  }
  legend(xy_legend[1], xy_legend[2], colnames(timemat),color)
}


################ INTERNAL FUNCTIONS ####################

ablineperso=function (a = NULL, b = NULL, h = NULL, v = NULL, reg = NULL, 
                      coef = NULL, untf = FALSE, x1 = NULL, x2 = NULL, y1 = NULL, 
                      y2 = NULL, ...) 
{
  if (!is.null(c(x1, x2, y1, y2))) {
    oldclip <- par("usr")
    if (is.null(x1)) 
      x1 <- oldclip[1]
    if (is.null(x2)) 
      x2 <- oldclip[2]
    if (is.null(y1)) 
      y1 <- oldclip[3]
    if (is.null(y2)) 
      y2 <- oldclip[4]
    clip(x1, x2, y1, y2)
    abline(h = oldclip[4] + 1)
    clip(x1, x2, y1, y2)
  }
  abline(a = a, b = b, h = h, v = v, reg = reg, coef = coef, 
         untf = untf, ...)
  if (!is.null(c(x1, x2, y1, y2))) 
    do.call("clip", as.list(oldclip))
}

# Return first common element to lists
ancestor <- function(lists){
  res = lists[[1]]
  for (i in lists[2:length(lists)]) {
    res = intersect(res, i)
  }
  return(res[1])
}

# Return list of all the ancestors of a node
ancestors <- function(node,gtreel){
  res=vector("numeric",0)
  ind<-node
  is_root=is.na(gtreel$nodes[ind,1])
  if (is.na(gtreel$nodes[node,1])){
    is_root=TRUE
  }
  while (is_root==FALSE){
    ind=gtreel$nodes[ind,1]
    res[length(res)+1]=ind
    is_root=is.na(gtreel$nodes[ind,1])
  }
  return(res)
}

# look for node of equal heights in the tree, if found that break the try by randomly changing the corresponding node heights
# write a set of trees if numrep>1
break_node_height_tie <- function(gtreel, tiplabel, numrep=100){
  inode_height = internal_node_height(gtreel)
  if ( sum( duplicated(inode_height[which(inode_height>0)]) )==0 ) {stop('no pair of internal nodes with equal heights')} # no internal nodes have identical heights
  node_id = sort(inode_height,decreasing=TRUE,index.return=TRUE)$ix # treat nodes from oldest to most recent
  treeset <- vector("list", numrep)
  class(treeset) <- "multiPhylo"
  for (nrep in 1:numrep){
    #     print(nrep)
    gtreelN = gtreel
    n=1
    while (inode_height[node_id[n]]>0) { # stop when reached the first tip
      #     x11(); plotree(gtreel,seq(1,15))
      m=n+1
      if (inode_height[node_id[n]]==inode_height[node_id[m]]) { # tie found need to break it
        #         print(paste("node:",node_id[n]))
        L_to_break=inode_height[node_id[n]]
        maxL=inode_height[node_id[n-1]]
        while (inode_height[node_id[m]]==L_to_break) {
          m=m+1
          minL=inode_height[node_id[m]]
        }
        #       print(c(node_id[n],maxL,minL))
        for (p in n:(m-1)){
          #         print(node_id[p])
          # sample a truncated normal between the height of previous and following nodes: Compute the CDF for the upper and lower limits of the interval and 
          # generate a uniform random numbers within that range of CDF values, then compute the inverse CDF
          # standard deviation is set as 5% of interval of height
          sdev=0.05*(maxL-minL)
          new_height = qnorm(runif(1, pnorm(minL, mean=L_to_break, sd=sdev), pnorm(maxL, mean=L_to_break, sd=sdev)), mean=L_to_break, sd=sdev)
          height_correction = new_height - inode_height[node_id[p]]
          #           print(height_correction)
          gtreelN$nodes[node_id[p],4] = gtreelN$nodes[node_id[p],4] + (height_correction)
          gtreelN$nodes[gtreelN$nodes[node_id[p],2:3],4] = gtreelN$nodes[gtreelN$nodes[node_id[p],2:3],4] - (height_correction)
        }
      }
      n=m
    }
    #     x11(); plotree(gtreelN,seq(1,15))
    treeset[[nrep]] = lconverttree(gtreelN)
  }
  return(treeset)
}

##check the structure of a gtreel phylogenetic tree
check_treel <- function(gtreel,locations=NA){
  error_found = 0
  if ( length(which(is.na(gtreel$nodes[,1])))!=1 ){
    print("error root")
    error_found = 1
  }
  for (n in 1:length(gtreel$nodes[,1])){
    pere = gtreel$nodes[n,1]
    filsg = gtreel$nodes[n,2]
    filsd = gtreel$nodes[n,3]
    if ( !is.na(pere) ){
      if ( gtreel$nodes[pere,2]!=n && gtreel$nodes[pere,3]!=n ){
        print("error 1")
        error_found = 1
      }
    }
    if (  (is.na(filsg) & !is.na(filsd))  |  (!is.na(filsg) & is.na(filsd))  ){
      print("error 2")
      error_found = 1
    }
    if ( !is.na(filsg) ){
      if ( gtreel$nodes[filsg,1]!=n ){
        print("error filsg")
        error_found = 1
      }
    }
    if ( !is.na(filsd) ){
      if ( gtreel$nodes[filsd,1]!=n ){
        print("error filsd")
        error_found = 1
      }
    }
    if (!is.na(locations[1])){
      if ( !is.na(filsg) & !is.na(filsd) ){
        if ( locations[filsg]!=locations[n] & locations[filsd]!=locations[n] ){
          print(paste("error locations, node",n))
          error_found = 1}
      }
    }
  }
  return(error_found)
}

##
### return the set of nodes below a given node
children <- function(node, gtreel){
  if (is.na(gtreel$nodes[node,2]) & is.na(gtreel$nodes[node,3])) {
    return(node)
  } else {
    return(c(node, children(gtreel$nodes[node,2], gtreel), children(gtreel$nodes[node,3], gtreel)))
  }
}

### return the set of tips below a given node
childrentip <- function(node, gtreel){
  if (is.na(gtreel$nodes[node,2]) & is.na(gtreel$nodes[node,3])) {
    return(node)
  } else {
    return(c(childrentip(gtreel$nodes[node,2], gtreel), childrentip(gtreel$nodes[node,3], gtreel)))
  }
}

## testing: return the time of coalescence for ind1 and ind2
coal_time <- function(coal_history,ind1,ind2){
  ev1 = NULL ## store all coalescent events for individual 1
  while(1){
    ev1_n = which(coal_history[,1]==ind1)
    if (length(ev1_n)==0){ev1_n = which(coal_history[,2]==ind1)}
    if (length(ev1_n)==0){break}
    else {
      ind1 = coal_history[ev1_n,3]
      ev1 = c(ev1,ev1_n)
    }
  }
  ev2 = NULL ## store all coalescent events for individual 2
  while(1){
    ev2_n = which(coal_history[,1]==ind2)
    if (length(ev2_n)==0){ev2_n = which(coal_history[,2]==ind2)}
    if (length(ev2_n)==0){break}
    else {
      ind2 = coal_history[ev2_n,3]
      ev2 = c(ev2,ev2_n)
    }
  }
  ## get the common coalescent events
  #print(ev1)
  #print(ev2)
  ev_co = intersect(ev1,ev2)[1]
  #print(ev_co)
  if (length(ev_co)==0){
    #print(0)
    #print('---')
    return(0)
  }else if (is.finite(ev_co)==FALSE){ ## is.finite(ev_co) test for NA
    #print(0)
    #print('---')
    return(0)
  }else{
    #print(coal_history[ev_co,4])
    #print('---')
    return(coal_history[ev_co,4])
  }
}

# Return the new gtreel #
conversion <- function(tree_phylo,corresp){
  treel_phylo=converttree(tree_phylo)
  label=(gsub("[a-z]","",tree_phylo$tip))
  label_int=sapply(label,as.integer,USE.NAMES=F)
  t = matrix(NA, length(treel_phylo[,1]),4)
  empty=c(1:length(t[,1]))
  full=which(!(is.na(t[,1])))
  void=(sum(!(is.na(empty)))==0)
  nodes=unique(t[which(!(is.na(t[,1])))])
  to_do=label_int
  root=vector('numeric',0)
  while(void==FALSE){
    for (i in (to_do)){
      label_c=which(corresp==i)
      ancestor_c=treel_phylo[label_c,1]
      child1_c=treel_phylo[label_c,2]
      child2_c=treel_phylo[label_c,3]
      t[i,1]=corresp[ancestor_c]
      t[i,2]=corresp[child1_c]
      t[i,3]=corresp[child2_c]
      t[i,4]=treel_phylo[label_c,4]
      
      if ((is.na(t[i,1])) & (!(is.na(t[i,2]))) ){
        root=i
      }
    }	
    full=which(!(is.na(t[,1])))
    if (length(root)==1){
      full[length(full)+1]=root			
    }
    empty[intersect(empty,full)]=NA
    nodes=unique(t[which(!(is.na(t[,1])))])
    to_do=setdiff(nodes,full)
    void=(sum(!(is.na(empty)))==0)
  }
  return(t)
}


### convert ape phylo tree structure (edges) to node index:parent|left child|right child|branch length to parent|
converttree <- function(tree_phylo){
  #gtreel = matrix(NA, (2*tree_phylo$Nnode + 1), 4)
  gtreel <- list( nodes=matrix(NA,(2*tree_phylo$Nnode+1),4), nodes.label=seq(1,(2*tree_phylo$Nnode+1)) )
  for (n in 1:(2*tree_phylo$Nnode)) {
    gtreel$nodes[tree_phylo$edge[n,2],1] = tree_phylo$edge[n,1] # parent
    if (is.na(gtreel$nodes[tree_phylo$edge[n,1],2])) { # has that edge been visited yet
      gtreel$nodes[tree_phylo$edge[n,1],2] = tree_phylo$edge[n,2] # child left
    } else {
      gtreel$nodes[tree_phylo$edge[n,1],3] = tree_phylo$edge[n,2] # child right
    }
    gtreel$nodes[tree_phylo$edge[n,2],4] = tree_phylo$edge.length[n]
    if (tree_phylo$edge[n,2]<=(tree_phylo$Nnode+1)) { gtreel$nodes.label[tree_phylo$edge[n,2]] = tree_phylo$tip.label[tree_phylo$edge[n,2]] } # tip
    if (!is.null(tree_phylo$node.label)) { gtreel$nodes.label[tree_phylo$edge[n,1]] = tree_phylo$node.label[tree_phylo$edge[n,1]-(tree_phylo$Nnode+1)] }
  }
  return(gtreel)
}

# Return tree without node labels 
correction <- function(tree_par){
  delete<-vector('numeric',0)
  ind=1
  t<-strsplit(tree_par,'')[[1]]
  for (i in 1:length(t)){
    if (t[i]==")"){
      n=(i+1)
      ponctuation=F
      while (ponctuation!=T){	
        delete[ind]=n
        ind=ind+1
        if (t[n+1]==";" | t[n+1]==":"){
          ponctuation=T
        }
        n=n+1					
      }
    }
  }
  for (i in delete){
    t[i]=NA
  }
  t<-(t[!is.na(t)])
  return(paste(unlist(t),collapse=''))
}

## Return the correspondence between gtreel and tree_phylo for tips and nodes #
correspondance <- function(tree_phylo,gtreel,change_topo=T){                    	
  treel_phylo=converttree(tree_phylo)
  label=(gsub("[a-z]","",tree_phylo$tip))
  label_int=sapply(label,as.integer,USE.NAMES=F)
  corresp=c(label_int,rep(NA,(length(gtreel[,1])-length(label_int))))
  if (change_topo==F){
    full=FALSE
    while (full==FALSE){
      for (i in (length(label_int)+1):length(gtreel[,1])){
        a=c(treel_phylo[i,2],treel_phylo[i,3])
        x=c(corresp[a[1]],corresp[a[2]])
        if ((!(is.na(x[1]))) & (!(is.na(x[2]))) ){
          corresp[i]=which(((gtreel[,2]==x[1]) | (gtreel[,3])==x[1]))					}
      }
      full=(sum(is.na(corresp))==0)
    }
  }else{
    all=c(1:length(gtreel[,1]))
    corresp[(length(label_int)+1):length(corresp)]=all[-label_int]
  }
  return(corresp)
}

### create a space by uniformly sampling coordinates
create_landscape <- function(space_size,space_dim,max_dist=10,usedegree=0){
  if(max_dist>20000)
  {
    stop("create_landscape(): Error maximum distance is greater than half equatorial circumference ")	
  }
  #If usedegree = 1 we use the Signed degrees format  
  if (space_dim==1) {
    space = runif(space_size-2,min=0,max=1)
    return(rbind(c(0,space,1)))
    #return(seq(0,1,length.out=space_size))
    # 	}else if (space_dim==2){
    # 		Lat = runif(space_size,min=0,max=max_dist)
    # 		Long = runif(space_size,min=0,max=max_dist)
    # 		mLat = max(Lat)-min(Lat)
    # 		mLong = max(Long)-min(Long)
    # 		if (mLat>mLong) {
    # 			Lat = (Lat-min(Lat))/mLat
    # 			Long = ( (Long-min(Long))/mLong ) * (mLong/mLat)
    # 		} else {
    # 			Lat = ( (Lat-min(Lat))/mLat ) * (mLat/mLong)
    # 			Long = (Long-min(Long))/mLong
    # 		}
    # 		return(rbind(Lat,Long))
    # 	}
  }else if (space_dim==2){
    if(usedegree==0){
      Lat = runif(space_size,min=0,max=max_dist)
      Long = runif(space_size,min=0,max=max_dist)
      return(rbind(Lat,Long))
    }
    else if(usedegree ==1){
      library(geosphere)
      #half_maxd=max_dist/2
      if(space_size==1)
      {
        return(rbind(lat1,long1))
      }
      else if(space_size>1)
      {
        Lat=c()
        Long=c()
        P1=c()
        long1=0
        lat1=runif(1,min = -66.56,max = 66.56) #artic & antartic circle
        P1[1]=long1
        P1[2]=lat1
        P2=c()
        #angle=runif(1,min = 0,max = 90)
        P2=destPoint(P1,45,max_dist*1000)
        long2=P2[1]
        lat2=P2[2]
        geopoints=matrix(nrow=space_size-2,ncol=2)
        n = 1
        while (n<=space_size-2){
          z = runif(1,min=-6371,max=6371)
          phi = runif(1,min=long1,max=long2)
          theta = asin(z/6371)
          lat = (180/pi)*theta#report our value on the y axis and runiform between both lat
          if ((lat1<lat2 && lat>=lat1 && lat<=lat2)||(lat2<lat1 && lat<=lat1 && lat>=lat2)){
            geopoints[n,]=c(phi,lat)
            n = n+1
          }
        }
        Long=c(Long,geopoints[,1])
        Lat=c(Lat,geopoints[,2])
        
        #Choose the starting point with marsaglia method
        # res1=c()
        # res2=c()
        # df=(matrix(rep(0,400), ncol = 4))
        # for(i in 1:100)
        # {
        #x1=runif(1,min = -1,max = 1)
        #  x2=runif(1,min = -1,max = 1)
        # while((x1^2+x2^2>=1))
        # {
        #   x1=runif(1,min = -1,max = 1)
        #   x2=runif(1,min = -1,max = 1)
        # }
        # x=2*x1 * sqrt(1-x1^2-x2^2)
        # y=2*x2 * sqrt(1-x1^2-x2^2)
        # z=1-2*(x1^2+x2^2)
        # 
        # r=sqrt(x^2+y^2+z^2)
        # teta=atan(y/x)
        # phi=acos(z/r)
        # long1=teta*180/pi
        # lat1=phi*180/pi
        # lat1=90-lat1
        # #long1=long1+90
        # #long1=long1*sample(c(1,-1),1)
        #   
        #   b=cos((lat1-90)*pi/180)
        #   c=tan(long1*pi/180)
        #   A=(2+2*(c^2)-2*(b^2)-2*(c^2)*(b^2))/(b^2)
        #   B=(4*(c^4)+8*(c^2)+3+(1/(b^2)))
        #   C=(4+4*(c^2))
        #   res1=c(res1,A*(x1^4)-B*(x1^2)+C)
        #   
        #   A=(-4*(b^2)*(x2^4))+(4*(x2^4))
        #   B=((x2^2)-(c^2))/2
        #   C=((c^4)*(b^2)-(c^4)-(b^2)-1)
        #   res2=c(res2,A+B+C)
        #   df[i,]=c(x1,x2,res1[i],res2[i])
        # }
        # 
        
        #Sample in a circle "http://gis.stackexchange.com/questions/25877/how-to-generate-random-locations-nearby-my-location"
        # for(i in 1:space_size-1){
        #   r=maxd*1000/2/111300
        #   u=runif(1, min = 0,max = 1)
        #   v=runif(1, min = 0,max = 1)
        #   w=r*sqrt(u)
        #   t=2*pi*v
        #   xeu=w*cos(t)
        #   yeu=w*sin(t)
        #   xprime=xeu/cos(lat1*pi/180)
        #   xtrouve=xprime+long1
        #   ytrouve=yeu+lat1
        #   Long=c(Long,xtrouve)
        #   Lat=c(Lat,ytrouve)
        # }
        
        
        #Use geosphere for all point
        # angle=45
        #  for(i in 1:(space_size-1))
        #  {
        #    dist_use=runif(1,min = 0,max = half_maxd)
        #    angle=runif(1, min = 0, max = 360)
        #    P2=destPoint(P1,angle,dist_use*1000)
        #    Long=c(Long,P2[1])
        #    Lat=c(Lat,P2[2])
        #  }
        
        #WRONG: Sample lat and long uniformely (does not take into account the curvature of the earth)
        # if(space_size>2)
        # {
        #   if(lat1<=lat2){
        #     Lat=runif(space_size-2,min=lat1,max=lat2)
        #   } else {
        #     Lat=runif(space_size-2,min=lat2,max=lat1)
        #   }
        #   
        #   if((long1<0 && long2>0) || (long2<0 && long1>0))
        #   {
        #     #calcule la somme des deux long 0 ou -180 et on regarde laquelle est la plus petite tha ton use pour les calcule
        #     if(long1<0)
        #     {
        #       neg=long1
        #       pos=long2
        #     }
        #     else 
        #     {
        #       neg=long2
        #       pos=long1
        #     }
        #     
        #     dist0=abs(abs(neg)+pos)
        #     dist180=abs((abs(neg)-180)+(180-pos))
        #     
        #     if(dist180<dist0)
        #     {
        #       if(long1<0)
        #       {
        #         dif1=abs(long1)-180
        #         dif2=180-long2
        #       }
        #       else
        #       {
        #         dif1=abs(long2)-180
        #         dif2=180-long1
        #       }
        #       pick=runif(space_size-2,min = dif1,max = dif2)
        #       for(i in 1:length(pick))
        #       {
        #         if(pick[i]<0){
        #           Long[i]=-180+abs(pick[i])
        #         }
        #         else{
        #           Long[i]=180-pick[i]
        #         }
        #       }
        #     }  
        #   }
        #   else if(long1<=long2){
        #     Long=runif(space_size-2,min=long1,max=long2)
        #   } else {
        #     Long=runif(space_size-2,min=long2,max=long1)
        #   }
        # }
        Lat=c(Lat,lat1,lat2)
        Long=c(Long,long1,long2)
        return(rbind(Lat,Long))
      }
    }
  }
}

# delete a migration in history, and all subsequent migrations that are affected
delete_mig_history <- function(history,mig_id_del) {
  #print(history)
  #print(mig_id_del)
  #print("---")
  loc_del = history[mig_id_del,2] # location to be discarded because not colonised anymore
  time2add = history[mig_id_del,4]
  if (mig_id_del<nrow(history)){
    for (mig_id in (mig_id_del:(nrow(history)-1))) {history[mig_id,]=history[mig_id+1,]} # delete this event from history
  }
  history = history[1:(nrow(history)-1),]
  history[mig_id_del,4] = history[mig_id_del,4]+time2add
  if (is.vector(history)) {return(history)} # need to keep at least one migration
  # need to delete migrations starting from discarded location, unless it has been colonised before or is location of root
  if ( (history[1,1]!=loc_del) && (loc_del %in% history[1:(mig_id_del-1),2])==FALSE){
    mig_id = mig_id_del
    while (mig_id<=nrow(history)){
      if (history[mig_id,1]==loc_del){
        history = delete_mig_history(history,mig_id)
        mig_id = mig_id-1 # need to stay on same index in case there has been deletion of migration(s) after
      }else if (history[mig_id,2]==loc_del) { # location being colonised again, no need to check anymore
        break
      }
      mig_id = mig_id+1
      if (is.vector(history)) {return(history)} # need to keep at least one migration
    }
  }
  return(history)
}

distkm <- function(lat1,lat2,long1,long2){ # get distance in km between coordinates in decimal (convert from degrees/minutes/seconds to decimal?)
    #d=sqrt(((long1-long2)*(long1-long2)) + ((lat1-lat2)*(lat1-lat2)))
    lat1r = (lat1/180) * pi # radian
    lat2r = (lat2/180) * pi
    long1r = (long1/180) * pi
    long2r = (long2/180) * pi
    R=6371
    x=(long2r-long1r)*cos(0.5*(lat2r+lat1r))
    y=lat2r-lat1r
    d=R*sqrt(x*x+y*y)
    if(d<=100)
    {
      return(d)
    }
    else
    {
      
      if(!identical(long1r,long2r) && abs(long2r-long1r)<1e-6){ #correspond to adistance lower than 1m
        d = acos( sin(lat1r)*sin(lat2r) + cos(lat1r)*cos(lat2r))*6378.137
        
      }else if(identical(long1,long2r) && identical(lat1r,lat2r))
      {
        d=0
      }
      else{
        #d=sqrt(((long1-long2)*(long1-long2)) + ((lat1-lat2)*(lat1-lat2)))
        d = acos( sin(lat1r)*sin(lat2r) + cos(lat1r)*cos(lat2r)*cos(long2r-long1r) ) * 6378.137
      }
    }
  return(d)
}


## compute effective sample size as N/(1+2*sum(autocorr))
effSaSize <- function(x){ ## http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introbayes_sect008.htm
  if (sum(x-x[1]) == 0) return(0) ## return 0 if all elements of x are equals
  autocor = acf(x,type="correlation",plot=FALSE,lag.max=length(x))[[1]]
  cutoff = which(autocor<(2*sd(x)))[1]
  if (is.na(cutoff)){
    return( NA )
  }else{
    return( length(x)/(1+2*sum(autocor[2:cutoff])) )
  }
}

# Return nucleotide frequencies for SeqGen 
freq_nucl <- function(upper_limit=0.5,lower_limit=0.1,save_freq=0){
  freq=rep(0.0,4)
  i=0
  while(i==0){
    frequences=rgamma(4,1)
    if(sum(frequences)!=0){
      freq=frequences/(sum(frequences))
    }
    i=1
    for (n in 1:length(freq)){
      if (freq[n]>upper_limit | freq[n]<lower_limit){i=0}
    }
  }
  if(save_freq==1){
    write.table(rbind(freq), file = "freq.txt", append = FALSE, quote = FALSE, sep = " ", na = "-", col.names=FALSE, row.names = FALSE,eol='')
  }
  return(freq)
}


### return the set of location of the tips below a node
ftips <- function(node, gtreel, loc_tips){
  if (is.na(gtreel$nodes[node,2]) & is.na(gtreel$nodes[node,3])) {
    return(loc_tips[node])
  } else {
    return(unique(c(ftips(gtreel$nodes[node,2], gtreel, loc_tips), ftips(gtreel$nodes[node,3], gtreel, loc_tips))))
  }
}

### return the set of location of the tips below a node WITH REPETITION
ftips_rep <- function(node, gtreel, loc_tips){
  if (is.na(gtreel[node,2]) & is.na(gtreel[node,3])) {
    return(loc_tips[node])
  } else {
    return(c(ftips_rep(gtreel[node,2], gtreel, loc_tips), ftips_rep(gtreel[node,3], gtreel, loc_tips)))
  }
}

## sample gene tree from species tree
geneal <- function(gtreel,numi,Npop=1){ ## numi is a vector of length of the number of species, giving the number of individuals for each species
  ## for example: numi=sample(1:10,20,replace=TRUE) to have between 1 and 10 individuals for each of the 20 species
  ## coalhist=geneal(gtreel,numi)
  ninternals = sum(numi)-1 ## number of leaves in the coalescence is sum(numi), -1 to have the number of internal nodes in the genealogy
  internal_id = (sum(numi)+1):(sum(numi)+ninternals)
  id_internal_id = 1
  coal_history = matrix(0,nrow=ninternals,ncol=4) ## store the history of coalescent events: node1|node2|internalnode|time
  ch_id = 1 ## store the index of coal_history
  uncoaled = vector("list",length(gtreel[,1])) ## store the list of un-coalesced individuals for each node in the species tree, if NULL then call indicoal()
  sroot = which(is.na(gtreel[,1])) ## root node in species tree
  snode_height = internal_node_height(gtreel) ## 
  ind_id = 1
  species = which(is.na(gtreel[,2]) & is.na(gtreel[,3])) ## tips of the species tree
  for (i in species){ ## get the list of individuals (therefore uncoalesced) for each tips in the species tree
    uncoaled[[i]] = ind_id:(ind_id + numi[i] - 1)
    ind_id = ind_id + numi[i]
  }
  #print(uncoaled)
  list_coal_history = recur_coal(gtreel,sroot,snode_height,coal_history,uncoaled,internal_id,id_internal_id,Npop) ## go through tree from tips to root, calling indicoal()
  coal_history = list_coal_history[[1]]
  #print(coal_history)
  uncoaled = list_coal_history[[2]]
  id_internal_id = list_coal_history[[3]]
  ## coalesce the remaining individuals at the root of the species tree
  if (length(uncoaled[[sroot]])>0){
    list_indicoal = indicoal(uncoaled[[sroot]],internal_id[id_internal_id:(id_internal_id+length(uncoaled[[sroot]])-2)],Npop)
    tmp_coal_history = list_indicoal[[1]]
    ## add up the height of the node to the coalescent events
    tmp_coal_history[,4] = tmp_coal_history[,4] + snode_height[sroot]
    if (coal_history[1,1]==0){
      id_coal_h = 1
    }else{
      id_coal_h = max(which(coal_history[,1]>0))+1 ## find the first coal event not recorded, "max(which(coal_history[,1]>0))" produces a warning when coal_history is empty (first instance) thus need the "if (coal_history[1,1]==0)"
    }
    coal_history[id_coal_h:(id_coal_h+length(tmp_coal_history[,3])-1),] = tmp_coal_history
  }
  ## reorder timewise the individual coalescent events history and clean up empty ones
  coal_history = coal_history[order(coal_history[which(coal_history[,1]>0),4]),]
  #print(uncoaled)
  #print(id_internal_id)
  return(list(uncoaled,coal_history))
}


## convert genealogy (history of coalescent events: |from|to|newnode|height|) to gtreel format (node index : |parent|left child|right child|branch length to parent|)
geneal_to_treel <- function(coal_history){
  gtreel = matrix(NA, length(unique(c(coal_history[,1],coal_history[,2],coal_history[,3]))), 4)
  for (n in 1:length(coal_history[,1])){
    p = coal_history[n,3]
    f1 = coal_history[n,1]
    f2 = coal_history[n,2]
    gtreel[f1,1] = p
    id_prev_f1 = which(coal_history[1:(n-1),3]==f1)
    if (length(id_prev_f1)>0){
      gtreel[f1,4] = coal_history[n,4]-coal_history[id_prev_f1,4]
    }else{
      gtreel[f1,4] = coal_history[n,4]
    }
    gtreel[f2,1] = p
    id_prev_f2 = which(coal_history[1:(n-1),3]==f2)
    if (length(id_prev_f2)>0){
      gtreel[f2,4] = coal_history[n,4]-coal_history[id_prev_f2,4]
    }else{
      gtreel[f2,4] = coal_history[n,4]
    }
    gtreel[p,2:3] = c(f1,f2)
  }
  return(gtreel)
}

### return the height of a node
height <- function(node, gtreel){
  if (is.na(gtreel$nodes[node,2]) & is.na(gtreel$nodes[node,3])) {
    return(0)
  } else {
    return(gtreel$nodes[gtreel$nodes[node,2],4] + height(gtreel$nodes[node,2], gtreel))
  }
}

## simulate coalescence of the individuals in list_ind, return a list of dated coalescent events
## if Npop==0, simulate a reverse Yule process with height of the tree equal to K
indicoal <- function(list_ind,internal_id,Npop,K=0){
  #print(paste(">",list_ind,"<"))
  if (length(internal_id)!=(length(list_ind)-1)) {stop("indicoal() error")}
  vect_ind = unlist(list_ind)
  tmp_uncoaled = vector("list",length(list_ind)-1) ## store the list of uncoalesced individuals at each step
  id_internal_id = 1
  tmp_coal_history = matrix(0,nrow=length(vect_ind)-1,ncol=4) ## store coalescent events history ind1|ind2|internal_id|time
  if (Npop==0){ ## Yule process, get factor for lambda of exp() so that height of the tree has mean of K
    lambda_factor = (1/K)*sum(1/(seq(1,(length(list_ind)-1),1)+1))
  }
  for (n in 1:(length(list_ind)-1)){
    tmp_coal_history[n,1:2] = sample(unique(vect_ind),2)
    if (n>1) {tlag=tmp_coal_history[n-1,4]}
    else {tlag=0}
    if (Npop>0){ ## Coalescence
      lambda = (length(unique(vect_ind))*(length(unique(vect_ind))-1))/(4*Npop) ## Felsenstein's book p456
    }else{ ## Yule process
      lambda = length(unique(vect_ind))*lambda_factor
    }
    tmp_coal_history[n,4] = rexp(1,rate=lambda) + tlag ## get the time of coalescence, prior to the last events (lag)
    #tmp_coal_history[n,4] = 0.0001+tlag ## DEBUG: so that all coalescent events are stored
    tmp_coal_history[n,3] = internal_id[id_internal_id] ## internal node they merged into
    vect_ind[which(vect_ind==tmp_coal_history[n,1] | vect_ind==tmp_coal_history[n,2])] = internal_id[id_internal_id] ## replace the coalesced individuals by the internal nodes they merged in
    tmp_uncoaled[[id_internal_id]] = unique(vect_ind)
    id_internal_id = id_internal_id + 1
  }
  #print(tmp_coal_history)
  return(list(tmp_coal_history,tmp_uncoaled))
}

### get all possible locations for internal nodes
internal_loc <- function(gtreel, loc_tips){
  locations = vector("list", length(gtreel$nodes[,1]))
  for (n in 1:length(gtreel$nodes[,1])) {
    if (is.na(gtreel$nodes[n,2]) & is.na(gtreel$nodes[n,3])) { ## tip
      locations[[n]] = loc_tips[n]
    } else {
      locations[[n]] = unique(c(ftips(gtreel$nodes[n,2], gtreel, loc_tips), ftips(gtreel$nodes[n,3], gtreel, loc_tips)))
    }
  }
  return(locations)
}

### get all possible locations for internal nodes WITH REPETITION
internal_loc_rep <- function(gtreel, loc_tips){
  locations = vector("list", length(gtreel[,1]))
  for (n in 1:length(gtreel[,1])) {
    if (is.na(gtreel[n,2]) & is.na(gtreel[n,3])) { ## tip
      locations[[n]] = loc_tips[n]
    } else {
      locations[[n]] = c(ftips_rep(gtreel[n,2], gtreel, loc_tips), ftips_rep(gtreel[n,3], gtreel, loc_tips))
    }
  }
  return(locations)
}

### get the node absolute heights of internal nodes from branch length to parent
internal_node_height <- function(gtreel){
  
  inode_height = vector("numeric", length(gtreel$nodes[,1]))
  
  for (n in 1:length(gtreel$nodes[,1])) {
    if (is.na(gtreel$nodes[n,2]) & is.na(gtreel$nodes[n,3])) { ## tip
      inode_height[n] = 0
    } else {
      inode_height[n] = gtreel$nodes[gtreel$nodes[n,2],4] + height(gtreel$nodes[n,2], gtreel)
    }
  }
  return(inode_height)
}

### convert louis tree structure to ape phylo tree structure (edges), assuming tips appear first in gtreel
lconverttree <- function(gtreel,alphab=0){
  tipid = which(is.na(gtreel$nodes[,2]))
  if (alphab==1){  tipid = as.integer(sort(as.character(tipid))) } #tips sorted alphabetically
  n = length(tipid)
  Nnode = n - 1
  nodeid = rep(0, Nnode)
  nodeid[1] = which(is.na(gtreel$nodes[,1])) # root
  edge = matrix(NA, nrow(gtreel$nodes)-1, 2)
  edgelength = vector('numeric', nrow(gtreel$nodes)-1)
  edgeid = 1
  for (i in 1:nrow(gtreel$nodes)){
    if (is.na(gtreel$nodes[i,1])) {i = i + 1 ; next} # root
    x = which(nodeid == gtreel$nodes[i,1])
    if (length(x) == 0){
      x = match(0, nodeid) # next free nodelabel
      nodeid[x] = gtreel$nodes[i,1]
    }
    edge[edgeid,1] = x + n
    z = which(tipid == i)
    if (length(z) == 0){ # is it a tip?
      y = which(nodeid == i)
      if (length(y) == 0){
        y = match(0, nodeid) # next free nodelabel
        nodeid[y] = i
      }
      edge[edgeid,2] = y + n
    }else{
      edge[edgeid,2] = z
    }
    edgelength[edgeid] = gtreel$nodes[i,4]
    edgeid = edgeid + 1
    i = i + 1
  }
  #tree_phylo = list(edge = edge, tip.label = paste("t",tipid,sep=""), node.label = paste("t", nodeid, sep=""), Nnode = Nnode, edge.length = edgelength)
  tree_phylo = list(edge = edge, tip.label = gtreel$nodes.label[tipid], node.label = gtreel$nodes.label[nodeid], Nnode = Nnode, edge.length = edgelength)
  class(tree_phylo) = "phylo"
  return(tree_phylo)
}

##
### convert louis tree structure to ape phylo tree structure (edges)
### optional tips sorted alphabetically
lconverttree2 <- function(gtreel,alphab=0){
  tiplabel = which(is.na(gtreel[,2]))
  if (alphab==1){	tiplabel = as.integer(sort(as.character(tiplabel))) } #tips sorted alphabetically
  n = length(tiplabel)
  Nnode = n - 1
  nodelabel = rep(0, Nnode)
  nodelabel[1] = which(is.na(gtreel[,1])) # root
  edge = matrix(NA, nrow(gtreel)-1, 2)
  edgelength = vector('numeric', nrow(gtreel)-1)
  edgeid = 1
  for (i in 1:nrow(gtreel)){
    if (is.na(gtreel[i,1])) {i = i + 1 ; next} # root
    x = which(nodelabel == gtreel[i,1])
    if (length(x) == 0){
      x = match(0, nodelabel) # next free nodelabel
      nodelabel[x] = gtreel[i,1]
    }
    edge[edgeid,1] = x + n
    z = which(tiplabel == i)
    if (length(z) == 0){ # is it a tip?
      y = which(nodelabel == i)
      if (length(y) == 0){
        y = match(0, nodelabel) # next free nodelabel
        nodelabel[y] = i
      }
      edge[edgeid,2] = y + n
    }else{
      edge[edgeid,2] = z
    }
    edgelength[edgeid] = gtreel[i,4]
    edgeid = edgeid + 1
    i = i + 1
  }
  tree_phylo = list(edge = edge, tip.label = paste("t",tiplabel,sep=""), node.label = paste("t", nodelabel, sep=""), Nnode = Nnode, edge.length = edgelength)
  class(tree_phylo) = "phylo"
  return(tree_phylo)
}

### convert louis tree structure to ape phylo tree structure (edges)
### WATCH OUT change the index of the root before doing the conversion
lconverttree_old <- function(gtreel,loc_sampled=0,internlab=0){
  ## test if the tips come first
  ## ...
  edge = matrix(NA, nrow(gtreel)-1, 2)
  elength = vector('numeric',nrow(gtreel)-1)
  tips = which(is.na(gtreel[,2]))
  numnodes = length(tips)-1
  ## root must be first internal node after tips
  if (which(is.na(gtreel[,1]))!=(max(tips)+1)){
    rootn = which(is.na(gtreel[,1]))
    new_rootn = max(tips)+1
    tmp = gtreel[rootn,]
    gtreel[rootn,] = gtreel[new_rootn,]
    gtreel[new_rootn,] = tmp
    ## swap the indexes
    gtreel[which(gtreel==new_rootn)] = 0
    gtreel[which(gtreel==rootn)] = new_rootn
    gtreel[which(gtreel==0)] = rootn
    if (loc_sampled[1]>0) {
      tmp = loc_sampled[new_rootn,]
      loc_sampled[new_rootn,] = loc_sampled[rootn,]
      loc_sampled[rootn,] = tmp
    }
  }
  id = 1
  for (n in 1:nrow(gtreel)) {
    if (n==which(is.na(gtreel[,1]))) {
      next
    } else {
      edge[id,1] = gtreel[n,1]
      edge[id,2] = n
      elength[id] = gtreel[n,4]
      id = id+1
    }
  }
  if (loc_sampled[1]>0) {## insert locations and individuals names at tips
    if (internlab==1){
      tree_phylo = list(edge=edge, edge.length=elength, tip.label=as.character(paste('L',loc_sampled[1:max(tips),1],'_i',tips,sep="")), Nnode=numnodes, node.label=as.character(paste('L',loc_sampled[new_rootn:(new_rootn+numnodes-1),1])))
    }else{## no internal labels
      tree_phylo = list(edge=edge, edge.length=elength, tip.label=as.character(paste('L',loc_sampled[1:max(tips),1],'_i',tips,sep="")), Nnode=numnodes, node.label=vector("character",numnodes))
    }
  }else{## no tip labels
    if (internlab==1){
      tree_phylo = list(edge=edge, edge.length=elength, tip.label=as.character(seq(1,max(tips),1)), Nnode=numnodes, node.label=as.character(seq(1,numnodes,1)))
    }else{## no internal labels
      tree_phylo = list(edge=edge, edge.length=elength, tip.label=as.character(seq(1,max(tips),1)), Nnode=numnodes, node.label=vector("character",numnodes))
    }
  }
  #print(tree_phylo)
  class(tree_phylo) = "phylo"
  return(tree_phylo)
}

##
### convert louis tree structure to ape phylo tree structure (edges), assuming tips appear first in gtreel
lconverttree_old2 <- function(gtreel){
  tiplabel = which(is.na(gtreel[,2]))
  if ( length(tiplabel) != ((nrow(gtreel)+1)/2) ) {
    print("lconverttree(): error gtreel format, uncorrect number of tips")
  }else if ( sum( tiplabel!=1:((length(gtreel[,1])+1)/2) ) ){
    print("lconverttree(): error gtreel format, tips do not appear first in gtreel")
  }
  n = length(tiplabel)
  Nnode = n - 1
  nodelabel = rep(0, Nnode)
  nodelabel[1] = which(is.na(gtreel[,1])) # root
  edge = matrix(NA, nrow(gtreel)-1, 2)
  edgelength = vector('numeric', nrow(gtreel)-1)
  edgeid = 1
  for (i in 1:nrow(gtreel)){
    if (is.na(gtreel[i,1])) {i = i + 1 ; next} # root
    edge[edgeid,] = c(gtreel[i,1],i)
    edgelength[edgeid] = gtreel[i,4]
    #nodelabel[edgeid] = gtreel[i,1]
    if (i>(n+1)) nodelabel[edgeid-n+1] = i # tips appear first in gtreel, followed by root, therefore first internal node is gtreel(n+1+1)
    edgeid = edgeid + 1
  }
  tree_phylo = list(edge = edge, tip.label = paste("t",tiplabel,sep=""), node.label = paste("t", nodelabel, sep=""), Nnode = Nnode, edge.length = edgelength)
  class(tree_phylo) = "phylo"
  return(tree_phylo)
}

## read MCMC output files and create data structure for plotting
# txtfile: data from simulations
# logfile: locations of nodes 
# csvfile: statistics on posterior
#txtfile="/home/louis/Documents/FRDF_steph_2010/Simulations/Retrieved/11Mar2012/phyloland_Sun11Mar2012_120243.9347_sim.txt"
#logfile="/home/louis/Documents/FRDF_steph_2010/Simulations/Retrieved/11Mar2012/phyloland_Sun11Mar2012_120243.9347_loc.log"
#csvfile="/home/louis/Documents/FRDF_steph_2010/Simulations/Retrieved/11Mar2012/phyloland_Sun11Mar2012_120243.9347.csv"
load_mcmc <- function(txtfile,logfile,csvfile){
  ## load simulations parameters
  con = file(txtfile, 'r')
  input = readLines(con, n=-1)
  close(con)
  space = rbind( as.numeric(unlist(strsplit(input[2]," "))), as.numeric(unlist(strsplit(input[3]," "))) )
  if (!require("ape",character.only=TRUE)){library(ape,character.only=TRUE)}
  tree_phylo = read.tree(text=input[5])
  #tree_phylo = reorder(tree_phylo, order = "cladewise")
  ## indexes have changed, need to update the tree so that indexes are read in node names ("tX")
  tree_phylo2 = tree_phylo
  tree_phylo2$edge = matrix(0,nrow(tree_phylo$edge),2)
  newidx = as.numeric(gsub("[a-z]","",tree_phylo$tip))
  for (n in 1:length(tree_phylo$tip)){ # first the tips, index equal reference in $edge
    m = which(tree_phylo$edge[,1]==n)
    tree_phylo2$edge[m,1] = newidx[n]
    m = which(tree_phylo$edge[,2]==n)
    tree_phylo2$edge[m,2] = newidx[n]
  }
  newidx = as.numeric(gsub("[a-z]","",tree_phylo$node))
  for (n in 1:length(tree_phylo$node)){ # second internal nodes, need to add the number of tips to get index in $edge
    m = which(tree_phylo$edge[,1]==(n+length(tree_phylo$tip)))
    tree_phylo2$edge[m,1] = newidx[n]
    m = which(tree_phylo$edge[,2]==(n+length(tree_phylo$tip)))
    tree_phylo2$edge[m,2] = newidx[n]
  }
  gtreel = converttree(tree_phylo2)
  if ( length(which(is.na(gtreel[,1])))>1 ){
    print('bug')
  }
  locations_sim = cbind( as.numeric(unlist(strsplit(input[7]," "))), 0, 0) # need three columns for compatibility
  ## load sampled node locations
  con = file(logfile, 'r')
  input = readLines(con, n=-1)
  close(con)
  input = strsplit(input[3:length(input)],"\t") # remove headers and first null line
  output_loc = do.call(rbind,lapply(input,as.numeric))
  sampledloc = output_loc[,2:length(output_loc[1,])] # remove first number which is the sample number
  ## load stat on posterior
  data = read.table(csvfile,sep = ",",header = FALSE,skip = 0,strip.white = TRUE)
  #stat_loc = sampled_mig(sampledloc,locations_sim,gtreel,space,1)
  if (length(which(is.na(gtreel[,1])))>1){print("error tree format");return()}
  height_dist = showlocat(sampledloc,gtreel,locations_sim,space)
  return(list(sampledloc,locations_sim,gtreel,space,data,height_dist))
}

## read MCMC output files from a directory
load_mcmc_dir <- function(dirnam){
  csvfiles = list.files(path=dirnam,pattern="phyloland[[:alnum:],_,.]*.csv",full.names=TRUE)
  tab_height_dist = matrix(0,0,2)
  for (n in 1:length(csvfiles)){
    csvfile = csvfiles[n]
    print(csvfile)
    txtfile = paste(sub(".csv","",csvfile),"_sim.txt",sep="")
    logfile = paste(sub(".csv","",csvfile),"_loc.log",sep="")
    if (!file.exists(txtfile) || !file.exists(logfile)){
      print(paste("output files not found for ",csvfile,sep=""))
    } else {
      tab_height_dist = rbind(tab_height_dist,load_mcmc(txtfile,logfile,csvfile)[[6]])
    }
  }
  x11();hist(tab_height_dist[which(tab_height_dist[,1]>0)])
  x11();plot(tab_height_dist[which(tab_height_dist[,1]>0),1],tab_height_dist[which(tab_height_dist[,1]>0),2])
  x11();cor(tab_height_dist[which(tab_height_dist[,1]>0),1],tab_height_dist[which(tab_height_dist[,1]>0),2])
}

# Return location of a common ancestor
loc_ancestor <- function(gtreel,tips,sampled_loc){
  anc=ancestor(lapply(X=tips,FUN=ancestors,gtreel=gtreel))
  return(as.integer(sampled_loc[anc]))
}

# improved list of objects
ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand
lsos <- function(..., n=10) {
    ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

##
## run MCMC for estimating variance and proba of colonising a previously occupied location
# freq: how often the estimated values are recorded
# Nstep_sigma: how many draws of sigma values at each iteration
# Nstep_lambda: how many draws of lambda at each iteration
# Nstep_loc: how many draws of internal locations at each iteration
# Nstep_genealogy: how many draws of genealogy at each iteration
# est_sigma: do not estimate sigma (0), use a uniform prior (1) or a dual prior (2)
mcmc_phyloland <- function(space, gtreel, simul_values, treelikelihood, est_sigma=c(1,1), est_lambda=1, est_tau=1, sample_geneal=0, sample_loc=0, 
                           plot_phylo=0, plot_MCMC=0, save_MCMC=0,
                           Nstep=1000, freq=1, Nstep_sigma=1, Nstep_lambda=1, Nstep_loc=1, Nstep_genealogy=1, Nstep_tau=1, pchange=0, ess_lim=100,
                           model=4, show_loc=0, filena_loc=NA, file_tracer=NA, file_tree=NA, dmethod="distkm"){
  
  if (!is.na(file_tracer)){
    write.table(t(c(0, rep(0, dim(space)[1]), 0, 0, 0)), file = file_tracer, row.names = FALSE, col.names = c("Sample", paste("Sigma", 1:dim(space)[1]), "Lambda", "Tau", "LogLikelihood"), sep = "\t", append = FALSE)
  }
  if (!is.na(file_tree)){
    #write("", file = file_tree, append = FALSE)
    file.create(file = file_tree)
  }
  if (!is.na(filena_loc)){
    file.create(file = filena_loc)
  }
  
  sigma_simul = simul_values[[1]]
  lambda_simul = simul_values[[2]]
  tau_simul = simul_values[[3]]
  locations_sim = simul_values[[4]]
  
  Ndim = dim(space)[1]
  tnumi = sum(is.na(gtreel$nodes[,2])) ## number of tips
  no_metropolis = 0 ## debug: to not use metropolis ratio (no data)
  
  #mat_Dists = space_dist(space)
  mat_Dists = space_dist(space, dmethod)
  
  if (sum(est_sigma)>0) sigma = rep(0,Ndim)
  print("Compute sigma upper limit(s)...")
  sig_lower_limit=c()
  sig_upper_limit=c()
  for (nd in 1:Ndim){
    sig_lower_limit=c(sig_lower_limit,sigma_limit_V2(dist =mean(mat_Dists[[nd]][mat_Dists[[nd]]>0]),proba_lim = 0.01 ))
    sig_upper_limit=c(sig_upper_limit,sigma_limit_V2(dist =mean(mat_Dists[[nd]][mat_Dists[[nd]]>0]),proba_lim = 0.99 ))
  }
  
  #sig_lower_limit = 0 # note that when new value of sigma proposed is very low then likelihood is -Inf then exp(metropolis+hastings)==0 then new value is rejected
  #sig_upper_limit = unlist(lapply(lapply(mat_Dists,max),sigma_limit,tries=1e3)) # limit on prior for sampling sigma during mcmc
  #sigma_threshold = unlist(lapply(lapply(mat_Dists,mean),sigma_limit)) # threshold below which limited dispersal is inferred
  sigma_threshold = sig_lower_limit
  if (sum(est_sigma==2)>0){
    sig_lam = rep(0,Ndim)
    sig_toplim = rep(0,Ndim)
    alpha = .5 ## probability for drawing prior on sigma equal to Uniform
  }
  for (nd in 1:Ndim){
    if (est_sigma[nd]>0){
      sigma[nd] = sig_upper_limit[nd]/2 
      print(paste("sigma upper limit",nd,":",format(round(sig_upper_limit[nd], 2), nsmall = 2)))
    }else{
      sigma[nd] = sigma_simul[nd]
    }
    if(est_sigma[nd]==2){
      sig_lam[nd] = -log(0.01)/sig_upper_limit[nd]
      sig_toplim[nd] = 1-exp(-sig_lam[nd]*sig_upper_limit[nd])
    }
  }
  
  lam_lower_limit = 1e-3
  lam_upper_limit = 1e3
  if (est_lambda==1){ ## start at random
    #lambda = runif(1,min=lam_lower_limit,max=lam_upper_limit)
    lambda = 1 ## start at 1, i.e. no competition
  }else{
    lambda = lambda_simul
  }
  
  tau_lower_limit = 1e-3
  tau_upper_limit = 1e3
  if (est_tau==1){ ## start at random
    #tau = runif(1,min=tau_lower_limit,max=tau_upper_limit)
    # start at average rate according to tree height and number of internal nodes
    tau = max(internal_node_height(gtreel)) / (tnumi-1)
  }else{
    tau = tau_simul
  }
  
  if (sample_loc==1 | is(sample_geneal,"multiPhylo")){ ## sample locations of all internal nodes, if sample_loc==0 the true internal nodes locations are used
    browser()
    if (is(sample_geneal,"multiPhylo")){
      ind = sample(1:length(sample_geneal),1)
      tree_phylo = sample_geneal[[ind]]
      gtreel = converttree(tree_phylo)
      #gtreelikelihood_curr = treelikelihood[[ind]]
    }else{
      tree_phylo = lconverttree(gtreel)
    }
    #replace real_loc by locations_sim (dans la fonction show_locat jesais pas d'ou vien le real_loc)
    dim_location=dim(locations_sim)[2]
    if(!is.null(dim_location))
    {
      locations_sim = locations_sim[,1]
    }
    locations = internal_loc(gtreel,locations_sim) ## get list of possible locations for all nodes, Uniform
    ## draw migrations uniformly
    locations_sampled = sim_history(gtreel,space,sig_upper_limit,1,tau_upper_limit/2,locations,0,model,1)[[1]]
    check_treel(gtreel,locations_sampled)
    if (show_loc==1){
      sort_space_id = sort(space,decreasing=TRUE,index.return=TRUE)$ix
      par(mfrow=c(2,1));
      xbar = barplot(hist(locations_sampled[,1],breaks=0:ncol(space),plot=FALSE)$counts[sort_space_id])
      title('true internal locations')
      histloc = matrix(nrow=(tnumi-1),ncol=Nstep) ## store location of internal nodes
      histdraw = vector("list",Nstep) ## store which internal node to be updated
    }
  }else{
    locations_sampled = locations_sim
    tree_phylo = lconverttree(gtreel)
  }
  
  ### acceptance rates computed on 100 values, first line values, second line acceptance
  acc_lambda = vector('numeric',100)
  id_acc_lambda = 1
  acc_sigma = matrix(0,Ndim,100)
  id_acc_sigma = vector('numeric',Ndim) + 1
  acc_tau = vector('numeric',100)
  id_acc_tau = 1
  acc_geneal = vector('numeric',100)
  id_acc_geneal = 1
  
  #### trees_phylo sampled
  trees_sampled = rmtree(1,2) # start by creating 1 random tree, 2 tips 
  list_trees = rmtree(1,2)
  list_trees = list_trees[-1]
  list_nb = vector("numeric",0)
  list_ind = vector("list",Nstep/freq)
  
  ### tuning parameters
  target_rate = 0.3 # acceptance rate to achieve
  init_tuningp = 1e2 ## start high and lower it progressively
  tuningp_lambda = init_tuningp
  tune_tp_lambda = 1 ## need to tune tuning parameter for lambda?
  tuningp_tau = init_tuningp
  tune_tp_tau = 1 ## need to tune tuning parameter for tau?
  tuningp_sigma = rep(0,Ndim) + init_tuningp
  tune_tp_sigma = rep(1,Ndim) ## need to tune tuning parameter for sigma?
  
  freq_n = 1
  freq_n1 = 1
  output_val = matrix(nrow=Nstep/freq,ncol=(4+Ndim))
  output_loc = matrix(nrow=Nstep/freq,ncol=(length(gtreel$nodes[,1])+1))
  
  ## compute likelihood of this particular phylogeny
  likeli_curr = sim_history(gtreel,space,sigma,lambda,tau,locations_sampled[,1],1,model)[[2]]
  likeli_new = likeli_curr # initialisation
  print(likeli_curr) # debug
  
  rate_acc = c(0,0,0)
  
  for (n in 1:Nstep) {
    
    ## sample internal location and/or genealogy
    
    if (is(sample_geneal,"multiPhylo")){
      for (n2 in 1:Nstep_genealogy) {
        ind = sample(1:length(sample_geneal),1)
        tree_phylo_new = sample_geneal[[ind]]
        gtreel_new = converttree(tree_phylo_new)
        gtreelikelihood_new = treelikelihood[[ind]]
        locations = internal_loc(gtreel_new,locations_sim[,1]) ## get list of possible locations for all nodes, uniform
        locations_sampled_new = sim_history(gtreel_new,space,sigma,lambda,tau,locations,0,model,1)[[1]]
        check_treel(gtreel,locations_sampled_new)
        likeli_new = sim_history(gtreel_new,space,sigma,lambda,tau,locations_sampled_new[,1],1,model)[[2]]
        metropolis = likeli_new + gtreelikelihood_new - likeli_curr - gtreelikelihood_curr
        accepted = 0
        if (runif(1,min=0,max=1)<exp(metropolis)){
          tree_phylo = tree_phylo_new
          gtreel = gtreel_new
          locations_sampled = locations_sampled_new
          likeli_curr = likeli_new
          gtreelikelihood_curr = gtreelikelihood_new
          accepted = 1
        }
        if (plot_phylo==1){
          gtreeape = lconverttree(gtreel)
          plot.phylo(gtreeape,show.node.label=TRUE)
          readline()
        }
      }
    } else if(sample_loc==1) { ## sample location of internal nodes
      for (n2 in 1:Nstep_loc) {
        locations_sampled_new = sim_history(gtreel,space,sigma,lambda,tau,locations,0,model,1)[[1]]
        check_treel(gtreel,locations_sampled_new)
        likeli_new = sim_history(gtreel,space,sigma,lambda,tau,locations_sampled_new[,1],1,model)[[2]]
        metropolis = likeli_new - likeli_curr
        if (is.na(metropolis)) {
          if (is.infinite(likeli_new)) {print("Infinite Likelihood 1")}
          else{print("metropolis error (sigma)")}
        } else {
          if (runif(1,min=0,max=1)<exp(metropolis)){
            locations_sampled = locations_sampled_new
            likeli_curr = likeli_new
          }
        }
      }
    }
    
    ## sample sigma (one per space dimension)
    
    for (nd in 1:Ndim){
      if (est_sigma[nd]>0){
        for (n2 in 1:Nstep_sigma) {
          if (est_sigma[nd]==1){
            sigma_new = sigma
            accepted = 0
            #print(sigma_new)
            sigma_new[nd] = sigma[nd] * exp(tuningp_sigma[nd]*(runif(1,min=0,max=1)-.5))
            #print(sigma_new)
            if (sigma_new[nd]>sig_lower_limit && sigma_new[nd]<sig_upper_limit[nd]) {
              likeli_new = sim_history(gtreel,space,sigma_new,lambda,tau,locations_sampled[,1],1,model)[[2]]
              metropolis = likeli_new - likeli_curr
              if (no_metropolis==1) {metropolis = 0}
              hastings = log(sigma_new[nd]) - log(sigma[nd])
              if (is.na(metropolis)) {
                if (is.infinite(likeli_new)) {print("Infinite Likelihood 2")}
                else{print("metropolis error (sigma)")}
              } else if (is.na(hastings)) {
                print("hastings error (sigma)")
              } else {
                if (runif(1,min=0,max=1)<(exp(metropolis+hastings))) {
                  sigma = sigma_new
                  likeli_curr = likeli_new
                  accepted = 1
                }
              }
            }
          } else if (est_sigma[nd]==2){
            sigma_new = sigma
            accepted = 0
            if (runif(1,min=0,max=1)<alpha){
              #sigma_new[nd] = 2.527
              sigma_new[nd] = sig_upper_limit[nd]
            }else{
              #sigma_new[nd] = -log(1-runif(1,min=0,max=0.99))/1.822386
              sigma_new[nd] = -log(1-runif(1,min=0,max=sig_toplim[nd]))/sig_lam[nd]
            }
            likeli_new = sim_history(gtreel,space,sigma_new,lambda,tau,locations_sampled[,1],1,model)[[2]]
            if (runif(1,min=0,max=1)<exp(likeli_new-likeli_curr)){
              sigma = sigma_new
              likeli_curr = likeli_new
              accepted = 1
            }
          }
          acc_sigma[nd,id_acc_sigma[nd]] = accepted
          id_acc_sigma[nd] = id_acc_sigma[nd] + 1
          if (id_acc_sigma[nd]==100){
            if (est_sigma[nd]==1 && tune_tp_sigma[nd]==1){ ## adjust tuning parameter
              accept_rate = sum(acc_sigma[nd,])/100
              if ( accept_rate<(target_rate*.5) ){
                tuningp_sigma[nd] = tuningp_sigma[nd] * 0.75
              } else if ( accept_rate<(target_rate*.9) ){
                tuningp_sigma[nd] = tuningp_sigma[nd] * 0.5
              } else if ( accept_rate<(target_rate*1.1) && accept_rate>(target_rate*.9) ){
                tune_tp_sigma[nd] = 0 ## stop tuning when reached approximately target_rate for the first time
                print(paste("tuning parameter sigma",nd,":",tuningp_sigma[nd]))
              } else if ( accept_rate>(target_rate*1.5) ){
                tuningp_sigma[nd] = tuningp_sigma[nd] * 1.5
              } else if ( accept_rate>(target_rate*1.1) ){
                tuningp_sigma[nd] = tuningp_sigma[nd] * 1.25
              }
            }
            id_acc_sigma[nd] = 1
          }
        }
        if (nd==2) rate_acc = rbind(rate_acc,c(sigma_new[nd],accepted,likeli_new)) # debug
      }
    }
    
    ## sample lambda
    if (est_lambda==1){
      for (n2 in 1:Nstep_lambda) {
        accepted = 0
        lambda_new = lambda * exp(tuningp_lambda*(runif(1,min=0,max=1)-.5))
        if (lambda_new>lam_lower_limit && lambda_new<lam_upper_limit) {
          likeli_new = sim_history(gtreel,space,sigma,lambda_new,tau,locations_sampled[,1],1,model)[[2]]
          metropolis = likeli_new - likeli_curr
          if (no_metropolis==1) {metropolis = 0}
          hastings = log(lambda_new) - log(lambda)
          if (is.na(metropolis)) {
            if (is.infinite(likeli_new)) {print("Infinite Likelihood 3")}
            else{print("metropolis error (lambda)")}
          } else if (is.na(hastings)) {
            print("hastings error (lambda)")
          } else {
            if (runif(1,min=0,max=1)<(exp(metropolis+hastings))) {
              lambda = lambda_new
              likeli_curr = likeli_new
              accepted = 1
            }
          }
        }
        acc_lambda[id_acc_lambda] = accepted
        id_acc_lambda = id_acc_lambda + 1
        if (id_acc_lambda==100){
          if (tune_tp_lambda==1){
            accept_rate = sum(acc_lambda)/100
            if ( accept_rate<(target_rate*.5) ){
              tuningp_lambda = tuningp_lambda * 0.75
            } else if ( accept_rate<(target_rate*.9) ){
              tuningp_lambda = tuningp_lambda * 0.5
            } else if ( accept_rate<(target_rate*1.1) && accept_rate>(target_rate*.9) ){
              tune_tp_lambda = 0 ## stop tuning when reached approximately 0.3 for the first time
              print(paste("tuning parameter lambda :",tuningp_lambda))
            } else if ( accept_rate>(target_rate*1.5) ){
              tuningp_lambda = tuningp_lambda * 1.5
            } else if ( accept_rate>(target_rate*1.1) ){
              tuningp_lambda = tuningp_lambda * 1.25
            }
            id_acc_lambda = 1
          }
        }
      }
    }
    
    ## sample tau
    if (est_tau==1){
      for (n2 in 1:Nstep_tau) {
        accepted = 0
        tau_new = tau * exp(tuningp_tau*(runif(1,min=0,max=1)-.5))
        if (tau_new>tau_lower_limit && tau_new<tau_upper_limit) {
          likeli_new = sim_history(gtreel,space,sigma,lambda,tau_new,locations_sampled[,1],1,model)[[2]]
          metropolis = likeli_new - likeli_curr
          if (no_metropolis==1) {metropolis = 0}
          hastings = log(tau_new) - log(tau)
          if (is.na(metropolis)) {
            if (is.infinite(likeli_new)) {print("Infinite Likelihood 4")}
            else{print("metropolis error (tau)")}
          } else if (is.na(hastings)) {
            print("hastings error (tau)")
          } else {
            if (runif(1,min=0,max=1)<(exp(metropolis+hastings))) {
              tau = tau_new
              likeli_curr = likeli_new
              accepted = 1
            }
          }
        }
        acc_tau[id_acc_tau] = accepted
        id_acc_tau = id_acc_tau + 1
        if (id_acc_tau==100){
          if (tune_tp_tau==1){
            accept_rate = sum(acc_tau)/100
            if ( accept_rate<(target_rate*.5) ){
              tuningp_tau = tuningp_tau * 0.5
            } else if ( accept_rate<(target_rate*.9) ){
              tuningp_tau = tuningp_tau * 0.75
            } else if ( accept_rate<(target_rate*1.1) && accept_rate>(target_rate*.9) ){
              tune_tp_tau = 0 ## stop tuning when reached approximately 0.3 for the first time
              print(paste("tuning parameter tau :",tuningp_tau))
            } else if ( accept_rate>(target_rate*1.5) ){
              tuningp_tau = tuningp_tau * 1.5
            } else if ( accept_rate>(target_rate*1.1) ){
              tuningp_tau = tuningp_tau * 1.25
            }
            id_acc_tau = 1
          }
        }
      }
    }
    
    ## record value
    if (n==(freq_n*freq)) {
      ## save tree, reorder before so that internal nodes are from oldest to most recent
      ordered = reorder_treel(gtreel,locations_sampled[,1])
      locations_ordered = ordered[[2]]
      tree_phylo_ordered = reorder_treep(tree_phylo)
      trees_sampled[[n/freq]] <- tree_phylo_ordered
      
      test <- test_tree(tree_phylo_ordered, list_trees, list_nb, list_ind, n ,freq)
      list_trees = test[[1]]
      list_nb = test[[2]]
      list_ind = test[[3]]
      
      if (!is.na(file_tree)){
        #print(tree_phylo_ordered$edge)
        #x11();plot(tree_phylo_ordered)
        #print(file_tree)
        write(write.tree(tree_phylo_ordered,file=""), file = file_tree, append = TRUE)
      }
      
      freq_n = freq_n + 1
      output_val[n/freq,1] = n
      
      ## save parameters
      for (idn in 1:Ndim){
        output_val[n/freq,(idn+1)] = sigma[idn]
      }
      
      output_val[n/freq,(idn+2)] = lambda
      output_val[n/freq,(idn+3)] = tau
      output_val[n/freq,(idn+4)] = likeli_curr
      stringo = sprintf("%i    ",output_val[n/freq,1])
      
      for (idn in 2:length(output_val[1,])){
        stringo = paste(stringo,sprintf("%.4f    ",output_val[n/freq,idn]))
      }
      print(stringo)
      
      ## save locations
      output_loc[n/freq,1] = n
      output_loc[n/freq,2:(length(gtreel$nodes[,1])+1)] = locations_ordered
      if (!is.na(filena_loc)){
        write.table(t(c(n,locations_ordered)), col.names = FALSE, row.names = FALSE, file = filena_loc, sep = "\t", eol = "\n", append = TRUE)
      }
      
      if (!is.na(file_tracer)){
        write.table(t(output_val[n/freq,]), file = file_tracer, row.names = FALSE, col.names = FALSE, sep = "\t", append = TRUE)
      }
      
      ## check effective sample size
      if (n==(10*freq_n1*freq)){
        freq_n1 = freq_n1 + 1
        ess_val = vector('numeric',(Ndim+2)) ## sigmas plus lambda and tau
        ness = 1
        for (idn in 2:(Ndim+3)){
          val = output_val[is.finite(output_val[,idn]),idn]
          ess_val[ness] = effSaSize(val)
          for (nd in 1:Ndim){
            if (est_sigma[nd] == 0){
              ess_val[nd] = ess_lim	
            }
          }
          if (est_lambda == 0){
            ess_val[3] = ess_lim	
          }
          if (est_tau == 0){
            ess_val[4] = ess_lim	
          }
          ness = ness+1
        }
        print(paste('ess:',toString(ess_val)))
        if (sum(ess_val<ess_lim)==0){
          L = sum(is.finite(output_val[,2]))
          output_val = output_val[1:L,]
          output_loc = output_loc[1:L,]
          break
        }
      }
    }
  }
  
  if (sample_loc==1 && show_loc==1){
    return(list(output_val,histloc,histdraw,tnumi,locations_sim,locations_sampled,space))
  }else{
    return(list(output_val, output_loc, trees_sampled, list_trees, list_nb, list_ind[1:length(list_trees)], rate_acc, sigma_threshold))
  }
}

# PLOT MCMC ESTIMATE vs TRUE, with .05% and .95% intervals
# system("ccato phyloland ~/Documents/phylogeo/test_package_louis/simul_Feb2013/jobs")
# data = read.table("~/Documents/phylogeo/test_package_louis/simul_Feb2013/jobs",sep = ",",header = FALSE,skip = 0,strip.white = TRUE)
mcmc_plot <- function(data,print_out=0){
  maxsigma = quantile(c(data[,28],data[,39]),probs=.99,names=FALSE)
  # SIGMA1
  x11()
  #postsig = (1-data[,19])*data[,28] + data[,19]*2.527 # stat*2.527 + (1-stat)*mean()
  #postsig05 = (1-data[,19])*data[,21] + data[,19]*2.527
  #postsig95 = (1-data[,19])*data[,25] + data[,19]*2.527
  
  # fitmeanest = loess(data[,30] ~ data[,2]) #31
  # fit05est = loess(data[,23] ~ data[,2]) #24
  # fit95est = loess(data[,27] ~ data[,2]) #28
  
  fitmeanest = loess(data[,76] ~ data[,74])
  fit05est = loess(data[,24] ~ data[,3])
  fit95est = loess(data[,28] ~ data[,3])
  simlam = seq(0,maxsigma,.5)
  preddmean <- predict(fitmeanest, simlam, se=TRUE)
  predd05 <- predict(fit05est, simlam, se=TRUE)
  predd95 <- predict(fit95est, simlam, se=TRUE)
  #plot(data[,30] ~ data[,2], ylab="estimation", xlab="simulation", xlim=c(0,maxsigma),ylim=c(0,maxsigma),main="sigma1",pch=20,cex=.5)
  plot(data[,30] ~ data[,2], ylab="estimation", xlab="simulation", xlim=c(0,maxsigma),ylim=c(0,maxsigma),main="sigma1")
  plot(data[,76] ~ data[,74], ylab="estimation", xlab="simulation", xlim=c(0,maxsigma),ylim=c(0,maxsigma),main="sigma1")
  points(simlam,simlam,type="l",col="red")
  lines(simlam,preddmean$fit, lty="solid", col="red", lwd=2)
  lines(simlam,predd05$fit, lty=2)
  lines(simlam,predd95$fit, lty=2)
  if (print_out==1){
    dev.print(file="./sigma1_estim.eps",height=10)
  }
  # SIGMA2
  x11()
  #postsig = (1-data[,30])*data[,39] + data[,30]*2.527 # stat*2.527 + (1-stat)*mean()
  #postsig05 = (1-data[,30])*data[,32] + data[,30]*2.527
  #postsig95 = (1-data[,30])*data[,36] + data[,30]*2.527
  fitmeanest = loess(data[,41] ~ data[,3]) #42
  fit05est = loess(data[,34] ~ data[,3]) #35
  fit95est = loess(data[,38] ~ data[,3]) #39
  
  # fitmeanest = loess(data[,42] ~ data[,4]) #42
  # fit05est = loess(data[,35] ~ data[,4]) #35
  # fit95est = loess(data[,39] ~ data[,4]) #39
  simlam = seq(0,130,.5)
  preddmean <- predict(fitmeanest, simlam, se=TRUE)
  predd05 <- predict(fit05est, simlam, se=TRUE)
  predd95 <- predict(fit95est, simlam, se=TRUE)
  plot(data[,41] ~ data[,3], ylab="estimation", xlab="simulation", xlim=c(0,maxsigma),ylim=c(0,maxsigma),main="sigma2")
  # plot(data[,42] ~ data[,4], ylab="estimation", xlab="simulation", xlim=c(0,maxsigma),ylim=c(0,maxsigma),main="sigma2")
  points(simlam,simlam,type="l",col="red")
  lines(simlam,preddmean$fit, lty="solid", col="red", lwd=2)
  lines(simlam,predd05$fit, lty=2)
  lines(simlam,predd95$fit, lty=2)
  if (print_out==1){
    dev.print(file="./sigma2_estim.eps",height=10)
  }
  # LAMBDA
  x11()
  fitmeanest = loess(data[,51] ~ data[,4])
  fit05est = loess(data[,44] ~ data[,4])
  fit95est = loess(data[,48] ~ data[,4])
  # fitmeanest = loess(data[,52] ~ data[,5])
  # fit05est = loess(data[,45] ~ data[,5])
  # fit95est = loess(data[,49] ~ data[,5])
  simlam = seq(0,10,.01)
  preddmean <- predict(fitmeanest, simlam, se=TRUE)
  predd05 <- predict(fit05est, simlam, se=TRUE)
  predd95 <- predict(fit95est, simlam, se=TRUE)
  plot(data[,51] ~ data[,4], ylab="estimation", xlab="simulation", xlim=c(0,10),ylim=c(0,10),main="lambda")
  # plot(data[,52] ~ data[,5], ylab="estimation", xlab="simulation", xlim=c(0,10),ylim=c(0,10),main="lambda")
  points(simlam,simlam,type="l",col="red")
  lines(simlam,preddmean$fit, lty="solid", col="red", lwd=2)
  lines(simlam,predd05$fit, lty=2)
  lines(simlam,predd95$fit, lty=2)
  if (print_out==1){
    dev.print(file="./lambda_estim.eps",height=10)
  }
  # LOG-LAMBDA
  x11()
  fitmeanest = loess(log(data[,51]) ~ log(data[,4]))
  fit05est = loess(log(data[,44]) ~ log(data[,4]))
  fit95est = loess(log(data[,48]) ~ log(data[,4]))
  # fitmeanest = loess(log(data[,52]) ~ log(data[,5]))
  # fit05est = loess(log(data[,45]) ~ log(data[,5]))
  # fit95est = loess(log(data[,49]) ~ log(data[,5]))
  simlam = seq(log(0.01),log(10),.01)
  preddmean <- predict(fitmeanest, simlam, se=TRUE)
  predd05 <- predict(fit05est, simlam, se=TRUE)
  predd95 <- predict(fit95est, simlam, se=TRUE)
  limits = c(log(0.01),log(10))
  plot(log(data[,51]) ~ log(data[,4]), ylab="estimation", xlab="simulation", xlim=limits, ylim=limits, main="log(lambda)", xaxt="n", yaxt="n")
  # plot(log(data[,52]) ~ log(data[,5]), ylab="estimation", xlab="simulation", xlim=limits, ylim=limits, main="log(lambda)", xaxt="n", yaxt="n")
  tick_seq = seq(limits[1], limits[2], length.out=10)
  #axis(1, at=tick_seq, labels=round(exp(tick_seq),2) )
  #axis(2, at=tick_seq, labels=round(exp(tick_seq),2) )
  axis(1, at=tick_seq, labels=round(tick_seq,2) )
  axis(2, at=tick_seq, labels=round(tick_seq,2) )
  points(simlam,simlam,type="l",col="red")
  lines(simlam,preddmean$fit, lty="solid", col="red", lwd=2)
  lines(simlam,predd05$fit, lty=2)
  lines(simlam,predd95$fit, lty=2)
  if (print_out==1){
    dev.print(file="./loglambda_estim.eps",height=10)
  }
  # TAU
  x11()
  fitmeanest = loess(data[,61] ~ data[,5])
  fit05est = loess(data[,54] ~ data[,5])
  fit95est = loess(data[,58] ~ data[,5])
  # fitmeanest = loess(data[,62] ~ data[,6])
  # fit05est = loess(data[,55] ~ data[,6])
  # fit95est = loess(data[,59] ~ data[,6])
  simlam = seq(0,10,.01)
  preddmean <- predict(fitmeanest, simlam, se=TRUE)
  predd05 <- predict(fit05est, simlam, se=TRUE)
  predd95 <- predict(fit95est, simlam, se=TRUE)
  plot(data[,61] ~ data[,5], ylab="estimation", xlab="simulation", xlim=c(0,10),ylim=c(0,10),main="tau")
  # plot(data[,62] ~ data[,6], ylab="estimation", xlab="simulation", xlim=c(0,10),ylim=c(0,10),main="tau")
  points(simlam,simlam,type="l",col="red")
  lines(simlam,preddmean$fit, lty="solid", col="red", lwd=2)
  lines(simlam,predd05$fit, lty=2)
  lines(simlam,predd95$fit, lty=2)
  if (print_out==1){
    dev.print(file="./tau_estim.eps",height=10)
  }
}

## convert history of migrations to one possible tree in format louis
## ignore the last migration that happens when the observation are done
mig2treel <- function(history,plot_tree=0){
  nmig = nrow(history)-1 ## keep last migration for tip branch length
  gtreel = list(nodes=NA,nodes.label=NA)
  gtreel$nodes = matrix(0,1+2*nmig,4)
  location = vector('numeric',1+2*nmig)
  locations_sampled = matrix(0,nrow=1+2*nmig,ncol=ncol(history)-3) ## location(destination) | distance(s)
  istip = vector('numeric',1+2*nmig)
  ## root
  location[1] = history[1,1]
  locations_sampled[1,] = c(history[1,1],rep(0,ncol(history)-4))
  gtreel$nodes[1,] = c(NA,2,3,NA)
  istip[1] = 1
  for (i in 1:nmig){
    possibletips = which(location==history[i,1] & istip==1)
    nodeid = possibletips[sample(length(possibletips),1)] ## uniform to choose which tip the migration starts from
    gtreel$nodes[nodeid,2:3] = c(2*i,2*i+1)
    #time_to_father = runif(1,0:1) #rexp(1,rate=1)
    wait_time = history[i,4]
    gtreel$nodes[nodeid,4] = gtreel$nodes[nodeid,4] + wait_time
    #gtreel[2*i,1:4] = c(nodeid,NA,NA,wait_time)
    #gtreel[2*i+1,1:4] = c(nodeid,NA,NA,wait_time)
    gtreel$nodes[2*i,1:3] = c(nodeid,NA,NA)
    gtreel$nodes[2*i+1,1:3] = c(nodeid,NA,NA)
    istip[nodeid] = 0
    gtreel$nodes[which(istip==1),4] = gtreel$nodes[which(istip==1),4] + wait_time ## increase node father heights of the existing tips
    istip[2*i] = 1
    istip[2*i+1] = 1
    staygo = sample(1:2,2,replace=FALSE) ## uniform to choose which node stays on its father loc
    location[2*i] = history[i,staygo[1]]
    location[2*i+1] = history[i,staygo[2]]
    locations_sampled[2*i,1:(length(history[i,])-3)] = c(history[i,staygo[1]],history[i,5:length(history[i,])])
    locations_sampled[2*i+1,1:(length(history[i,])-3)] = c(history[i,staygo[2]],history[i,5:length(history[i,])])
    #nodeid = nodeid + 1
    #print(istip)
  }
  gtreel$nodes.label=location
  ## extra step to define the length of the tip branches, this is given by the last migration
  gtreel$nodes[which(istip==1),4] = gtreel$nodes[which(istip==1),4] + history[nmig+1,4]
  #print(istip)
  if (plot_tree==1){
    print(history)
    print(rbind(seq(1,length(location)),location))
    print(gtreel$nodes)
    treepar = treel2par(gtreel)
    #print(locations_sampled)
    #print(treepar)
    plot.phylo(read.tree(file="",text=treepar),show.tip.label=TRUE,show.node.label=TRUE)
    print(treepar)
  }
  # reformat gtreel for consistency
  #gtreelout = list(nodes=NA,nodes.label=NA)
  #gtreelout$nodes = gtreel
  return(list(gtreel ,location,locations_sampled))
}

### return a migration event
### occupied: >0=yes, 0=not
### if departure[1]>0 then only consider migrations from these locations
### if destination[1]>0 then only consider migrations to these locations
### if mig[1]>0 then only consider this particular migration out of all possible (must be vector of 3 elements with third element being the time to wait for this migration to happen)
migration<-function(space,occupied,sigma,lambda,tau,departure=0,destination=0,mig=0){
  #print(c(occupied,sigma,lambda,tau))
  if (length(sigma)!=dim(space)[1]) stop('migration(): error 1')
  if (sum(occupied)==0) stop('migration(): error 2')
  
  ## consider all events simultaneously, build a density matrix, normalise it and then sample an event in it, occupied locations can be re-colonised
  if (departure[1]>0) {
    departures = departure
  } else {## occupied locations
    departures = which(occupied>0)
  }
  if (destination[1]>0) {
    destinations = destination
  } else {## all locations
    destinations = 1:length(occupied)
  }
  
  if (dim(space)[1]==1) {## 1D
    mat_dists = abs(outer(space[departures],space[destinations],"-")) # SHOULD BE PUBLIC
    #print(mat_dists)
    log_densitymat = dnorm(mat_dists,mean=0,sd=sqrt(sigma),log=TRUE) - dnorm(0,mean=0,sd=sqrt(sigma),log=TRUE) # normalise by the density value for staying on the same location
    #print(log_densitymat)
    #print(occupied)
    #print(destinations)
    #print(departures)
  } else if (dim(space)[1]==2) {## 2D
    mat_distsL = abs(outer(space[1,departures],space[1,destinations],"-")) # SHOULD BE PUBLIC
    mat_distsl = abs(outer(space[2,departures],space[2,destinations],"-")) # SHOULD BE PUBLIC
    log_densitymat = (dnorm(mat_distsL,mean=0,sd=sqrt(sigma[1]),log=TRUE) - dnorm(0,mean=0,sd=sqrt(sigma[1]),log=TRUE)) + 
                     (dnorm(mat_distsl,mean=0,sd=sqrt(sigma[2]),log=TRUE) - dnorm(0,mean=0,sd=sqrt(sigma[2]),log=TRUE)) -
                     log(length(occupied))
    #log_densitymat = dnorm(mat_distsL,mean=0,sd=sigma[1],log=TRUE) + dnorm(mat_distsl,mean=0,sd=sigma[2],log=TRUE) #20Feb2013
    #densitymat = densitymat - ( log((pnorm(1,mean=0,sd=sigma[1])-pnorm(0,mean=0,sd=sigma[1]))) + log((pnorm(1,mean=0,sd=sigma[2])-pnorm(0,mean=0,sd=sigma[2]))) ) # not needed anymore since normalised by dnorm(0,...)
  }
  if (lambda==0){ # impossible to colonise an occupied location
    densitymat = exp(log_densitymat)
    densitymat[,which(occupied==1)] = 0
  }else{
    mat_occ = t( matrix(log(lambda)*(occupied[destinations]>0), nrow=length(destinations), ncol=length(departures)) )
    log_densitymat = log_densitymat + mat_occ
    densitymat = exp(log_densitymat)
  }
  densitymat[is.infinite(densitymat)] = 0 ## forbid -Inf values (impossible moves)
  if (sum(densitymat)==0) { ## no migration possible
    return(c(0,0,0,Inf,rep(0,dim(space)[1])))
  }
  #exp_param = rowSums(densitymat) * occupied[departures] * tau # factorise by the number of individuals occupying each location
  R = densitymat * tau
  Rs = 0
  for (n in 1:length(departures)){
    for (m in 1:length(destinations)){
      Rs = Rs + R[n,m]*(occupied[destinations[m]])
    }
  }
  p = rowSums(R)/Rs # proba to select edge to branch
  #q = matrix(0,length(departures),length(destinations))
  #for (n in 1:length(departures)){
  #  for (m in 1:length(destinations)){
  #    q[n,m] = R[n,m]/sum(R[n,]) # proba to select location to migrate to
  #  }
  #}
#print(exp_param)
  if (mig[1]>0) { ## find the migration event: from mig[1] to mig[2] and get probability of it
    lstart = which(departures==mig[1])[1]
    lgoto = which(destinations==mig[2])[1]
    ## taking branch length into account: proba_mig * proba_branch_length * proba_select_origin
    ## (the minimum of several exponential dist is distributed exp too, and the proba of the index is known, see wikipedia)
    #proba_event = densitymat[lstart,lgoto] * dexp(mig[3],sum(exp_param)) # drawing in n exponential and taking minimum is same as drawing in exponential which parameter is the sum of n parameters
    #proba_event = R[lstart,lgoto] * dexp(mig[3],Rs)
    if (sum(occupied)==1){ # this is the root, first migration event
      proba_event = R[lstart,lgoto] / sum(R[lstart,])
    }else{
      proba_event = R[lstart,lgoto] * exp(-Rs*mig[3])
    }
    wait_time = 0
  }else{ ## sample one migration event
    #     if (length(departures)>1){
    #   		wait_times = sapply(exp_param,function(x) rexp(1,x))
    #   		lstart = sort(wait_times,index.return=TRUE)$ix[1]
    #   	}else{ ## only one departure possible
    #   		wait_times = rexp(1,exp_param)
    #   		lstart = 1
    #   	}
    #     probavec = densitymat[lstart,]/sum(densitymat[lstart,])
    # #print(densitymat)
    # #print(lstart)
    # #print(probavec)
    # 		lgoto = sample(1:length(densitymat[lstart,]),1,prob=probavec)
    # 		proba_event = probavec[lgoto]
    # 		wait_time = wait_times[lstart]
    # sampling n exponential is = to sample one exponential which parameter is the sum
    wait_time = rexp(1,Rs)
    lstart = sample(1:length(p),1,prob=p)
    q = R[lstart,]/sum(R[lstart,]) # proba to select location to migrate to
    lgoto = sample(1:length(q),1,prob=q)
    #proba_event = densitymat[lstart,lgoto] * dexp(wait_time,sum(exp_param))
    #proba_event = p[lstart]*q[lgoto] ##20Feb2013
    proba_event = 0 #DEBUGGING return the likelihood of the migration event:
    if (sum(occupied)==1){ # this is the root, first migration event
      proba_event = R[lstart,lgoto] / sum(R[lstart,])
    }else{
      proba_event = R[lstart,lgoto] * exp(-Rs*wait_time)
    }
  }
  mig_event = vector('numeric',4+dim(space)[1])
  mig_event[1:4] = c( departures[lstart], destinations[lgoto], proba_event, wait_time )
  for (n in 1:dim(space)[1]) {
    mig_event[4+n] = abs(space[n,destinations[lgoto]]-space[n,departures[lstart]])
  }
  return(mig_event)
}

# wrapper function to call migration in C
### return a migration event
### occupied: >0=yes, 0=not
### if departure[1]>0 then only consider migrations from these locations in mode 2
### if destination[1]>0 then only consider migrations to these locations in mode 2
### if mig[1]>0 then only consider this particular migration out of all possible in mode 2 (must be vector of 3 elements with third element being the time to wait for this migration to happen)
new_migration <- function(space, occupied, sigma, lambda, tau, departure, destination, mig){
  migr <- function(space, occupied, sigma, lambda, tau, mig, mig_event, space_dim, space_size, length_mig){
    .C("migC", space, occupied, sigma, lambda, tau, mig, mig_event, space_dim, space_size, length_mig)}
  
  space_dimC = as.integer(nrow(space))
  space_sizeC = as.integer(ncol(space))
  
  spaceC = as.double(as.vector(t(space)))
  occupiedC = as.integer(occupied)
  sigmaC  = as.double(sigma)
  lambdaC = as.double(lambda)
  tauC = as.double(tau)
  migC = as.double(mig)
  
  #departure = as.integer(departure)
  #length_departure = as.integer(length(departure))
  
  #destination = as.integer(destination)
  #length_destination = as.integer(length(destination))
  
  mig_eventC = as.double(vector("numeric", (4+space_dimC)))
  length_migC = as.integer(length(migC))
  
  if (length(sigma) != space_dimC) stop('new_migration error 1')
  mig_event = migr(spaceC, occupiedC, sigmaC, lambdaC, tauC, migC, mig_eventC, space_dimC, space_sizeC, length_migC)[[7]]
  #mig_event = migration(space, occupied, sigma, lambda, tau, departure, destination, mig) ## DEBUG by using R version
  
  return(mig_event)
}

### print memory usage per variable
object.sizes <- function(){
  return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name) object.size(get(object.name))))))
}


output_function = function(stat_res,space_size,maxdist,ntip,sigma_simul,lambda_simul,tau_simul,est_sigma,est_lambda,est_tau,sample_geneal,sample_loc,Nstep,freq,ess_lim,sig_upper_limit,sig_lower_limit,space,gtreel,mean_matDist1,mean_matDist2,indice){
  # header for output log file
  headr = paste('num location',sep=",")
  headr = paste(headr,'num tip',sep=",")
  for (ndim in 1:2){
    headr = paste(headr,',sigma',ndim,sep="")
  }
  headr = paste(headr,'lambda','tau','sim mol seq','est_sigma1','est_sigma2','est_lambda','est_tau','sample_geneal','sample_loc','Nstep','freq','Nstep_sigma','Nstep_lambda','Nstep_tau','Nstep_loc','pchange', 'ess_lim',sep=",")
  for (ndim in 1:2){
    headr = paste(headr,paste('sigma',ndim,' stat',sep=""),'min','quant05','quant25','quant50','quant75','quant95','max','mode','mean','last',sep=",")
  }
  headr = paste(headr,'lambda_min','quant05','quant25','quant50','quant75','quant95','max_lambda','mode_lambda','mean_lambda','last_lambda',sep=",")
  headr = paste(headr,'tau_min','quant05','quant25','quant50','quant75','quant95','max_tau','mode_tau','mean_tau','last_tau',sep=",")
  for (ndim in 1:2){
    headr = paste(headr,paste('sim root loc coord ',ndim,sep=""),sep=",")
  }
  headr = paste(headr,',Nstep_genealogy',sep="")
  for (ndim in 1:2){
    headr = paste(headr,paste('sigma upper limit ',ndim,sep=""),sep=",")
  }
  for (ndim in 1:2){
    headr = paste(headr,paste('sigma lower limit ',ndim,sep=""),sep=",")
  }
  headr = paste(headr,paste('Max Dist ',sep=""),sep=",")
  for (ndim in 1:2){
    headr = paste(headr,paste('mean distance ',ndim,sep=""),sep=",")
  }
  id_filena = paste(format(Sys.time(),"%a%d%b%Y_%H%M%S"),sample(0:1000,1),sep="_")
  filena = paste('./phyloland_',indice,'_',id_filena,'.csv',sep='')
  cat(paste(headr,'\n',sep=""),file=filena,sep="",append=TRUE)
  filena_sim = paste('./phyloland_g_',id_filena,'_sim.txt',sep='')
  fileNWK=paste('tree_',id_filena,'.nwk',sep='')#     Nom du fichier contenant l'arbre au format parenthese
  fileNEX=paste('tree_seq_',id_filena,'.nex',sep='')# Nom du fichier contenant les sequences au format Nexus
  fileXML1=paste('file1_',id_filena,'.xml',sep='')#   Nom du fichier en entree de Beast au fornat xml
  fileXML2=paste('file2_',id_filena,'.xml',sep='')#   Nom du fichier en entree de Beast au fornat xml
  filena_loc = paste('./phyloland_g_',id_filena,'_loc1.log',sep='')
  #filena_acc_rate = paste('./phyloland_g_',id_filena,'_acc_rate1.txt',sep='')
  file_tracer = paste("tracer1_",id_filena,".log",sep="")
  file_tree = paste("file_tree1_",id_filena,".nwk",sep="")
  
  cat(c(space_size,ntip,round(sigma_simul,3),round(lambda_simul,3),round(tau_simul,3),0,est_sigma,est_lambda,est_tau,as.numeric(is(sample_geneal,"multiPhylo")),
        sample_loc,Nstep,freq,Nstep_sigma=1,Nstep_lambda=1,Nstep_tau=1,Nstep_loc=1,pchange=0,ess_lim,unlist(stat_res),round(space[,gtreel$nodes.label[which(is.na(gtreel$nodes[,1]))]],3),
        Nstep_genealogy=1,round(sig_upper_limit,3),round(sig_lower_limit,3),round(maxd,3),round(mean_matDist1,3),round(mean_matDist2,3),'\n'),file=filena,sep=",",append=TRUE)
  
  cat("space coordinates",'\n',sep="",file=filena_sim)
  cat("simulated tree",'\n',sep="",file=filena_sim,append=TRUE)
  gtreel$nodes[,4]=round(gtreel$nodes[,4],3)
  cat(treel2par(gtreel),'\n',sep=" ",file=filena_sim,append=TRUE)
  cat("simulated tree node locations",'\n',sep="",file=filena_sim,append=TRUE)
  cat(gtreel$nodes.label,'\n',sep=" ",file=filena_sim,append=TRUE)
  cat("simulated tree format louis",'\n',sep="",file=filena_sim,append=TRUE)
}

## permute tip locations of a tree
permute_loc <- function(gtreel,tiplabels1){
  tiplabels = tiplabels1[sample(length(tiplabels1))]
  locations = vector("list",length(gtreel[,1]))
  for (n in 1:length(gtreel[,1])) {
    if (is.na(gtreel[n,2]) & is.na(gtreel[n,3])) { ## tip
      #locations[[n]] = n
      locations[[n]] = as.integer(gsub("t","",tiplabels[n]))
    } else {
      locations[[n]] = c(ftips(gtreel[n,2],gtreel,tiplabels),ftips(gtreel[n,3],gtreel,tiplabels))
    }
  }
  return(locations)
}

## convert sequences from PHYLIP format to FASTA format
phylip_2_fasta <- function(filena){
  x = read.dna(filena)
  write.dna(x, paste(filena,".fasta",sep=""), format = "fasta")
}

# plot estimates from simulations, with spread between quantile estimates 5%-95%
plot_estim <- function(data,tp,mp,plot_mig=0){
  # tp: true value index for a given parameter
  # mp: median index for estimated values assuming quantile05 index is mp-2 and quantile95 index is mp+2
  #x11();par(mfrow=c(2,1));plot(data[,tp],data[,mp]);plot(data[,tp],abs(data[,mp+2]-data[,mp-2]));
  x11();
  intervls = hist(data[,tp],plot=FALSE,40)$breaks
  mean_spread = vector('numeric',length(intervls)-1)
  for (n in 2:length(intervls)){
    idx = which(data[,tp]>intervls[n-1] & data[,tp]<=intervls[n])
    if (length(idx)==0){
      mean_spread[n] = 0
    }else{
      mean_spread[n] = mean(abs(data[idx,mp+2]-data[idx,mp-2]))
    }
  }
  plot(data[,tp],data[,mp],xlab="",ylab="",xlim=c(min(data[,tp],intervls),max(data[,tp],intervls)),ylim=c(min(data[,mp],mean_spread),max(data[,mp],mean_spread)));
  #limites = par("usr")
  #print(limites)
  #par(new=TRUE);
  #plot(intervls,mean_spread,"l",xlab="",ylab="",xlim=c(limites[1],limites[2]),ylim=c(limites[3],limites[4]))#,axes=FALSE
  points(intervls,mean_spread,"l")
  #lines(mean_spread)
  #limites = par("usr")
  #axis(4,c(0,.5,1,1.5,2))
  #axis(1)
  #print(limites)
  if (plot_mig==1){
    # stats on root location
    sprintf("correlation root mean coord1: %.2f",cor(data[,61],data[,65]))
    sprintf("correlation root mean coord2: %.2f",cor(data[,62],data[,66]))
    sprintf("mean distance root coord. (mode): %.2f [%.2f,%.2f]",mean( sqrt((data[,61]-data[,63])^2+(data[,62]-data[,64])^2) ),
            quantile(sqrt((data[,61]-data[,63])^2+(data[,62]-data[,64])^2),.05,names=FALSE),
            quantile(sqrt((data[,61]-data[,63])^2+(data[,62]-data[,64])^2),.75,names=FALSE))
    sprintf("mean distance root coord. (mean): %.2f [%.2f,%.2f]",mean( sqrt((data[,61]-data[,65])^2+(data[,62]-data[,66])^2) ),
            quantile(sqrt((data[,61]-data[,65])^2+(data[,62]-data[,66])^2),.05,names=FALSE),
            quantile(sqrt((data[,61]-data[,65])^2+(data[,62]-data[,66])^2),.75,names=FALSE))
    ## get statistics on migration retrieval
    sensitiv = data[,68]/(data[,68]+data[,69]) # sim1sam1/(sim1sam1+sim1sam0)
    #data2 = data[which(data[,2]<0.5 & data[,3]<0.5),]
    #sensitiv = data[which(data[,2]<.5 & data[,3]<.5),68]/(data[which(data[,2]<.5 & data[,3]<.5),68]+data[which(data[,2]<.5 & data[,3]<.5),69])
    x11(); hist(sensitiv,xlim=c(0,1))
    sprintf("sensitivity: %.2f [%.2f,%.2f]",mean(sensitiv),quantile(sensitiv,.05,names=FALSE),quantile(sensitiv,.75,names=FALSE))
    specific = data[,71]/(data[,71]+data[,70]) # sim0sam0/(sim0sam0+sim0sam1)
    #specific = data[which(data[,2]<.5 & data[,3]<.5),71]/(data[which(data[,2]<.5 & data[,3]<.5),71]+data[which(data[,2]<.5 & data[,3]<.5),70])
    x11(); hist(specific,xlim=c(0,1))
    sprintf("specificity: %.2f [%.2f,%.2f]",mean(specific),quantile(specific,.05,names=FALSE),quantile(specific,.75,names=FALSE))
  }
}

plot_prior_posterior <- function(data,sig_lower_limit,sig_upper_limit,sigma_simul){
  for(i in 1:length(sig_lower_limit))
  {
    sigma=data[,i+1]
    # Get the density estimate
    #Calcule value pour echelle
    #hist(sigma,breaks = 10,xlim = c(0,2500000),ylim = c(0, 500))
    dens=density(sigma,from = sig_lower_limit[i], to = sig_upper_limit[i])
    # Plot y-values scaled by number of observations against x values
    plot(dens$x,dens$y,type="l",xlab="Sigma Value",ylab="Density",col = "blue")
    title(main = paste("Prior & Posterior distribution for Sigma ", i))
    abline(v=sigma_simul[i],col = 'red') 
    ablineperso(h = 1/(sig_upper_limit[i]-sig_lower_limit[i]), x1 = sig_lower_limit[i],x2 = sig_upper_limit[i],col = 'green')
    legend("topright", legend=c("Prior", "Posterior","True sigma"), col=c("green", "blue","red"), lty=1:2, cex=0.8, title="Distribution", text.font=3,bty = "n")
  }
}

#Function from the script with the KM modif
plot_sigma_limit <- function(sigma_min,sigma_max,dis_min,dis_max){
  xdim=seq(0,dis_max*1.1,length.out=100)
  plot(xdim,dnorm(xdim,0,sqrt(sigma_min))/(dnorm(0,0,sqrt(sigma_min))),type="l",lwd=2,col="blue",xlab="distance",ylab="probability density function")
  points(xdim,dnorm(xdim,0,sqrt(sigma_max))/(dnorm(0,0,sqrt(sigma_max))),type="l",lwd=2,col="red")
  #text(c(dis_min,dis_min), c(dnorm(dis_min*1.1,0,.3)/(pnorm(1,0,.3)-pnorm(0,0,.3)),dnorm(dis_min*1.1,0,sig_lower_limit)/(pnorm(1,0,sig_lower_limit)-pnorm(0,0,sig_lower_limit))),
  #     labels = c("0.3",sig_lower_limit))
  abline(v=dis_min, col="green")
  abline(v=dis_max, col="green")
  abline(h = 0.1)
  abline(h = 0.9)
}

# plot tree with location names for internal nodes and tips
# "location" contains the index of location from the oldest up to the most recent internal nodes, therefore needs to make sure tree_phylo is ordered the same way
plot_trees_tips<-function(tree_phylo, location, space, show_index=0){
  # make sure the internal nodes indexes are in correct order with location
  ordered = reorder_treel(converttree(tree_phylo),location)
  #location = ordered[[2]]
  tree_phylo = reorder_treep(tree_phylo)
  if (length(colnames(space))>0){ # locations of space are named
    nodes = colnames(space)
  }
  #   if (show_index==1){
  #     tips = paste(tree_phylo$tip.label, colnames(location)[as.numeric(sub("t","",tree_phylo$tip.label))], sep="_")
  #     tree_phylo$node.label = paste(as.numeric(location[(length(tree_phylo$tip.label)+1):(length(location))]), nodes[ as.numeric(location[(length(tree_phylo$tip.label)+1):(length(location))]) ], sep="_")
  #   }else{
  #     tips = colnames(location)[as.numeric(sub("t","",tree_phylo$tip.label))]
  #     tree_phylo$node.label = nodes[ as.numeric(location[(length(tree_phylo$tip.label)+1):(length(location))]) ]
  #   }
  tree_phylo$node.label = colnames(space)[ as.numeric(location[(Ntip(tree_phylo)+1):(2*Ntip(tree_phylo)-1)]) ]
  #  tree_phylo$tip.label = tips
  tree_par = write.tree(tree_phylo)
  plot.phylo(read.tree(file="", text=tree_par), show.tip.label=TRUE, show.node.label=TRUE, cex=.7, label.offset=.001)
  #plot.phylo(read.tree(file="", text=tree_par), show.tip.label=TRUE, show.node.label=TRUE, cex=.7) # Ape's label.offset option seems unstable??
  print(tree_par)	
}

### plot tree in format louis using ape package (or ade4)
plotree <- function(gtreel,location){
  require(gsubfn)
  if (missing("location")) location=paste("t",seq(1,nrow(gtreel$nodes)),sep="")
  treepar = treel2par(gtreel)
  ## replace node id by their location id (needs library(gsubfn))
  treepar = gsubfn('t([0-9]*):',function(x) paste(as.numeric(x),'-L',location[as.numeric(x)],':',sep=''),treepar, engine = "R")
  treepar = gsubfn('t([0-9]*);',function(x) paste(as.numeric(x),'-L',location[as.numeric(x)],';',sep=''),treepar, engine = "R") ## special case for the root
  plot.phylo(read.tree(file="",text=treepar),show.tip.label=TRUE,show.node.label=TRUE)
  print(treepar)
}

## recursively go through the tree
recur_coal <- function(gtreel,n,snode_height,coal_history,uncoaled,internal_id,id_internal_id,Npop){
  #print(paste("n=",n,"---"))
  if(length(uncoaled[[gtreel[n,2]]])==0){
    list_coal_history = recur_coal(gtreel,gtreel[n,2],snode_height,coal_history,uncoaled,internal_id,id_internal_id,Npop)
    coal_history = list_coal_history[[1]]
    uncoaled = list_coal_history[[2]]
    id_internal_id = list_coal_history[[3]]
  }
  if(length(uncoaled[[gtreel[n,3]]])==0){
    list_coal_history = recur_coal(gtreel,gtreel[n,3],snode_height,coal_history,uncoaled,internal_id,id_internal_id,Npop)
    coal_history = list_coal_history[[1]]
    uncoaled = list_coal_history[[2]]
    id_internal_id = list_coal_history[[3]]
  }
  ## left child
  if (length(uncoaled[[gtreel[n,2]]])==1){
    uncoaled[[n]] = uncoaled[[gtreel[n,2]]]
  } else {
    #print(tree[n,2])
    list_indicoal = indicoal(uncoaled[[gtreel[n,2]]],internal_id[id_internal_id:(id_internal_id+length(uncoaled[[gtreel[n,2]]])-2)],Npop)
    tmp_coal_history = list_indicoal[[1]]
    ## add up the height of the node to the coalescent events
    tmp_coal_history[,4] = tmp_coal_history[,4] + snode_height[gtreel[n,2]]
    tmp_uncoaled = list_indicoal[[2]]
    min_t = snode_height[n]
    coal_events_to_keep_ids = which(tmp_coal_history[,4]<min_t)
    #print(coal_events_to_keep_ids)
    #print(uncoaled)
    #print(uncoaled[[tree[n,2]]])
    #print(min_t)
    #print(n)
    #print(tmp_coal_history)
    if (length(coal_events_to_keep_ids)>0){
      tmp_coal_history = tmp_coal_history[coal_events_to_keep_ids,] ## select the coalescent events that happened before the age of the node
      if (length(coal_events_to_keep_ids)==1){tmp_coal_history = t(as.matrix(tmp_coal_history))}
      #print(tmp_coal_history)
      #print(">")
      #print(tmp_coal_history[,3])
      #print(">")
      if (coal_history[1,1]==0){
        id_coal_h = 1
      }else{
        id_coal_h = max(which(coal_history[,1]>0))+1 ## find the first coal event not recorded, "max(which(coal_history[,1]>0))" produces a warning when coal_history is empty (first instance) thus need the "if (coal_history[1,1]==0)"
      }
      coal_history[id_coal_h:(id_coal_h+length(tmp_coal_history[,3])-1),] = tmp_coal_history
      uncoaled[[n]] = tmp_uncoaled[[max(coal_events_to_keep_ids)]]
      id_internal_id = which(internal_id == tmp_coal_history[length(tmp_coal_history[,3]),3]) + 1 ## last internal node id used +1
    } else {
      uncoaled[[n]] = uncoaled[[gtreel[n,2]]]
    }
  }
  #print(uncoaled)
  #print(coal_history)
  #print(paste('---'))
  ## right child
  if (length(uncoaled[[gtreel[n,3]]])==1){
    uncoaled[[n]] = c(uncoaled[[n]],uncoaled[[gtreel[n,3]]])
  } else {
    #print(internal_id)
    #print(id_internal_id:(id_internal_id+length(uncoaled[[tree[n,3]]])-2))
    #print(tree[n,3])
    list_indicoal = indicoal(uncoaled[[gtreel[n,3]]],internal_id[id_internal_id:(id_internal_id+length(uncoaled[[gtreel[n,3]]])-2)],Npop)
    tmp_coal_history = list_indicoal[[1]]
    ## add up the height of the node to the coalescent events
    tmp_coal_history[,4] = tmp_coal_history[,4] + snode_height[gtreel[n,3]]
    tmp_uncoaled = list_indicoal[[2]]
    min_t = snode_height[n]
    coal_events_to_keep_ids = which(tmp_coal_history[,4]<min_t)
    #print(coal_events_to_keep_ids)
    #print(tmp_coal_history)
    if (length(coal_events_to_keep_ids)>0){
      tmp_coal_history = tmp_coal_history[coal_events_to_keep_ids,] ## select the coalescent events that happened before the age of the node
      if (length(coal_events_to_keep_ids)==1){tmp_coal_history = t(as.matrix(tmp_coal_history))}
      #print(tmp_coal_history)
      #print(">")
      #print(tmp_coal_history[,3])
      #print(">")
      if (coal_history[1,1]==0){
        id_coal_h = 1
      }else{
        id_coal_h = max(which(coal_history[,1]>0))+1 ## find the first coal event not recorded, max(which(coal_history[,1]>0)) produces a warning when coal_history is empty (first instance) thus need the "if (coal_history[1,1]==0)"
      }
      coal_history[id_coal_h:(id_coal_h+length(tmp_coal_history[,3])-1),] = tmp_coal_history
      uncoaled[[n]] = c(uncoaled[[n]],tmp_uncoaled[[max(coal_events_to_keep_ids)]])
      id_internal_id = which(internal_id == tmp_coal_history[length(tmp_coal_history[,3]),3]) + 1 ## last internal node id used +1
    } else {
      uncoaled[[n]] = c(uncoaled[[n]],uncoaled[[gtreel[n,3]]])
    }
  }
  #print(tmp_coal_history)
  #print(uncoaled)
  #print(paste('---'))
  #print(coal_history)
  return(list(coal_history,uncoaled,id_internal_id))
}



# read locations file (output interface)
read_loc <- function(fileLOC){
  sampled_loc <- read.table(file=fileLOC, header=FALSE, sep="\t", quote='"', comment.char="", check.names=FALSE)
  sampled_loc = sampled_loc[,-1]
  return(sampled_loc)
}

# read space file : a text file with 3 columns (separated by tabs). The first one with the tips names (the same as in \code{fileTrees}), the second one with the location latitudes and the last one with the location longitudes (both in decimal degrees). No header.
# return the tips locations, the space with latitudes and longitudes, the space [0,1]
read_space <- function(fileDATA, tips){
  sampled_loc <- as.matrix(read.table(fileDATA, header = FALSE, sep="\t"))
  if (length(sampled_loc)<(3*length(readLines(fileDATA)))){ # parsing failed, try sep=" " 
    sampled_loc <- as.matrix(read.table(fileDATA, header = FALSE, sep=" "))
  }
  tips_file = as.character(sampled_loc[,1])	
  if (sum(sort(tips_file) == sort(tips)) != length(tips)){
    print(sort(tips_file)[which(sort(tips_file) != sort(tips))])
    stop("in interface : tips names in fileDATA not in fileTREES",  call.=FALSE)	
  }
  coord_tips = matrix(c(as.double(sampled_loc[,2]),as.double(sampled_loc[,3])),ncol = 2)
  loc_tips <- vector("numeric",dim(coord_tips)[1])
  if (sum(duplicated(coord_tips[,1]) & duplicated(coord_tips[,2]))>0) { # remove duplicated locations from space
    space = t(coord_tips[-which(duplicated(coord_tips[,1]) & duplicated(coord_tips[,2])),])
  } else {
    space = t(coord_tips)
  }
  for (i in 1:length(loc_tips)){
    loc_tips[which(tips==tips_file[i])]  =  which(space[1, ] == coord_tips[i, 1] & space[2, ] == coord_tips[i, 2])
  }
  space_dec = space
  # normalise space in 0-1 in each dimension
  for ( i in 1:dim(space)[1]){
    L = space[i,]
    if (sum(L < 0) != 0){
      abso = abs(min(L))
      for (j in 1:length(L)){
        L[j] = L[j] + abso	
      }
    }
    mL = max(L) - min(L)
    space[i,] = (L - min(L))/mL	
  }	
  return(list(loc_tips,space,space_dec))	
}

# Read gtreel format 
read_treel <- function(gtreel){
  tr=as.matrix((read.table(gtreel,na.strings="-")))
  return(tr)
}

##read probabilities
read_treelikelihood <- function(fileTREES){
  conn=file(fileTREES,open="r")
  linn=readLines(conn) #read lines
  linn=grep("tree ",linn,value=T)#subset lines contain "tree"
  ##grep("tree +$",y,value=F) 
  B=length(linn)
  posterior=numeric(B) #create variable to hold posterior probabilities
  for (i in 1:length(linn)){
    if (length(linn)>0){
      #posterior[i]=unlist(strsplit( (unlist(strsplit(linn[i], split=c('posterior=',"]"), fixed=TRUE))[2]) ,split=']',fixed=T))[1]
      posterior[i]=unlist(strsplit( (unlist(strsplit(linn[i], split=c('posterior=',"]"), fixed=TRUE))[2]) ,split=',',fixed=T))[1]
    } ##substring
  }
  close(conn) 
  posterior<- grep("-[0-9]+\\.[0-9]+$",posterior,value=T) #get rid of results not containing -.[0-9]
  posterior2=numeric(B)  #create another variable to hold posterior probabilities if NO posterior probabilities in beast file. 
  for (i in 1:length(linn)){
    if (length(posterior)==0){
      posterior2[i]=1/length(linn) #posterior=1/Number of trees
    }
    else {
      posterior2=posterior #if posterior probabilities are provided from beast output
      treelikelihood=posterior2
    }
  }
  return(as.numeric(treelikelihood))
}


# reOrder tree format louis, so that oldest internal node first after the root
# tips must appear first, then root and finally internal nodes (case when gtreel is output by converttree())
reorder_treel <- function(gtreel,locations=0){
  gtreelO <- list( nodes=matrix(0,nrow=nrow(gtreel$nodes),ncol=ncol(gtreel$nodes)), nodes.label=rep(NA,length(gtreel$nodes.label)) )
  locationsO = locations*0 # intialise to 0
  ntips = ( length(gtreel$nodes[,1])+1 )/2
  nnodes = ntips-1
  if (!is.na(gtreel$nodes[ntips+1,1])){stop('reorder_treel() error: root not following tips')}
  gtreelO$nodes[1:ntips,2:3] = cbind(rep(NA,ntips),rep(NA,ntips)) # tips first (no children)
  gtreelO$nodes[ntips+1,c(1,4)] = NA # then root (no parent)
  gtreelO$nodes[1:ntips,4] = gtreel$nodes[1:ntips,4] # height for tips stay the same
  ## get the absolute height of each internal node
  inode_height = internal_node_height(gtreel)
  ## need to treat the nodes from the most recent to the oldest 
  node_id = sort(inode_height,decreasing=FALSE,index.return=TRUE)$ix
  rootn = which(is.na(gtreel$nodes[,1]))
  if (node_id[length(node_id)]!=rootn){stop('reorder_treel() error: root')} #wesh
  nh = cbind(node_id,c(node_id[1:ntips],nrow(gtreel$nodes):(rootn+1),rootn))
  locationsO[nh[,2]] = locations[nh[,1]]
  for (n in (ntips+1):nrow(gtreel$nodes)) {
    if (n==(ntips+1)) gtreelO$nodes.label[n] = gtreel$nodes.label[rootn]
    if(gtreel$nodes[nh[n,1],2]<=ntips){
      filsg = gtreel$nodes[nh[n,1],2]
    }else{
      filsg = nh[which(nh[,1]==gtreel$nodes[nh[n,1],2]),2]
    }
    if(gtreel$nodes[nh[n,1],3]<=ntips){
      filsd = gtreel$nodes[nh[n,1],3]
    }else{
      filsd = nh[which(nh[,1]==gtreel$nodes[nh[n,1],3]),2]
    }
    gtreelO$nodes[nh[n,2],2:3] = c(filsg,filsd)
    gtreelO$nodes[filsg,1] = nh[n,2]
    gtreelO$nodes[filsd,1] = nh[n,2]
    gtreelO$nodes[filsg,4] = gtreel$nodes[gtreel$nodes[nh[n,1],2],4]
    gtreelO$nodes[filsd,4] = gtreel$nodes[gtreel$nodes[nh[n,1],3],4]
    gtreelO$nodes.label[filsg] = gtreel$nodes.label[gtreel$nodes[nh[n,1],2]]
    gtreelO$nodes.label[filsd] = gtreel$nodes.label[gtreel$nodes[nh[n,1],3]]
  }
  return(list(gtreelO,locationsO))
}

# reOrder tree format Ape, so that oldest internal node appear first
# in Ape tree format, tips of the tree are always from 1 to number of tips of tree_phylo and internal nodes are indexed from Ntip(tree_phylo)+1
reorder_treep <- function(tree_phylo){
  tree_phyloO = tree_phylo
  ntips = Ntip(tree_phylo)
  edges = tree_phylo$edge
  edgesO = edges #rep(0,2*ntips-1)
  ## get the depth of each node
  nod_depth = node.depth.edgelength(tree_phylo)
  ## need to treat the internal nodes from the most recent to the oldest 
  nod_depth = nod_depth[(ntips+1):length(nod_depth)]
  node_id = sort(nod_depth,decreasing=TRUE,index.return=TRUE)$ix
  node_id = c(1:ntips,node_id+ntips)
  nh = cbind(node_id,c(node_id[1:ntips],(2*ntips-1):(ntips+1)))
  for (n in 1:length(nh[,1])){
    edgesO[which(edges==nh[n,1])] = nh[n,2]
  }
  tree_phyloO$edge = edgesO
  return(tree_phyloO)
}

## simulate migration in space and return corresponding genealogy
sim_phyloland <- function(N=10,Mean_H=1,pop_size=1,ni=10,space_size=10,sigma_simul=.4,lambda_simul=.5,tau_simul=1,use_molseq=0,space_dim=1,model=1,usedegree=0){
  if (model!=4){ ## create tree first
    ## get a tree using a Yule process then simulate a genealogy
    #tree_treeshape = rtreeshape(n=1, tip.number=N, model="yule")[[1]] # tree format "treeshape" (package:apTreeshape)
    #tree_phylo = as.phylo.treeshape(tree_treeshape) # tree format "phylo" (package:ape)
    #gtreel = converttree(tree_phylo) # tree format "louis": for each node index : |parent|left child|right child|branch length to parent|
    coal_species_hist = indicoal(seq(1,N),seq(N+1,2*N-1),0,Mean_H)
    gtreel = geneal_to_treel(coal_species_hist[[1]])
    #treep = write.tree(tree_phylo,file="",append=FALSE,digits=10,tree.names=FALSE) # parenthesis format
    #treep = gsub('tip_','t',treep) # tip names are t*
    ## simulate genalogy using multispecies coalescent
    ## PARAMETER: number of individuals for each species
    #numi = sample(1:ni,N,replace=TRUE)
    numi = rep(2,N)
    coalhist = geneal(gtreel,numi,pop_size)
    gtreel = geneal_to_treel(coalhist[[2]]) ## genealogy in format "louis"
  }else{ ## tree will be generated from migration events
    #gtreel = NA
    gtreel = list(nodes=NA,nodes.label=NA)
    
  }
  ## generate a space
  ## PARAMETER: number of geographic locations
  
  space = create_landscape(space_size,space_dim)
  
  #print(numi)
  #space = create_landscape(sum(numi),space_dim)
  #space = t(as.matrix(space))
  #print(space)
  ## simulate migrations from root to tips
  ## PARAMETER: sigma, lambda, tau
  simulated_history = sim_history(gtreel,space,sigma_simul,lambda_simul,tau_simul,0,0,model,0,N)
  
  locations_sampled = simulated_history[[1]]
  
  if (model==4){
    gtreel = simulated_history[[3]]
  }
  # ## get the list of possible locations for internal nodes including root
  #loctip = locations_sampled[1:length(which(is.na(gtreel[,2]) & is.na(gtreel[,3])))]
  #locations = internal_loc_rep(gtreel,loctip)
  #locations = vector("list",length(gtreel[,1]))
  #locations[which(is.na(gtreel[,2]) & is.na(gtreel[,3]))] = locations_sampled[which(is.na(gtreel[,2]) & is.na(gtreel[,3]))]
  #locations[which(!is.na(gtreel[,2]) & !is.na(gtreel[,3]))] = list(1:ncol(space))
  if (use_molseq==1){
    ## generate molecular sequences on the genealogy
    write.tree(lconverttree(gtreel,locations_sampled),file="./test1.tree")
    #system(paste("seq-gen -on -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -l40 < ./test1.tree > ./test1.dat",sep="")) ## HKY
    #system(paste("seq-gen -op -mHKY -l50 < ./test1.tree > ./test1.dat",sep="")) ## JC69
    #sim_seq = read.dna("./test1.dat")
    # ## need to convert sequences from PHYLIP to FASTA format (using {ape})
    #phylip_2_fasta("test1.dat")
    # ## multiple align sequences
    #system(paste("clustalw2 -INFILE=./test1.dat.fasta -OUTFILE=./test1.dat.fasta.nex -ALIGN -OUTPUT=NEXUS"))
    ## generate trees
    write_beauti(sim_seq,numi,"./test1.xml")
    system("beast -overwrite ./test1.xml")
    trees = read.nexus("./test1.trees") # list of trees of class "phylo"
  }
  if (use_molseq==0){
    # 		## delete unoccupied locations
    # print(space)
    # 		space = space[sort(unique(locations_sampled[,1]))] ## watch out, need to keep the smae order (hence sort())
    # 		space = (space-min(space))/(max(space)-min(space)) ## scale space between 0 and 1
    # 		space = t(as.matrix(space)) ## convert back to matrix for code consistency
    # print(space)
    # 		## replace index with new ones
    # 		unoccupied = vector('numeric',sum(numi))+1
    # 		unoccupied[unique(locations_sampled[,1])] = 0
    # print(unoccupied)
    # print(locations_sampled)
    # 		for (j in 1:max(locations_sampled[,1])){
    # 			locations_sampled[locations_sampled[,1]==j,1] = j-(sum(unoccupied[1:(j-1)]))
    # 		}
    # print(locations_sampled)
    return(list(space,gtreel,locations_sampled))
  }else{
    return()
  }
}

## get some stats on migration and locations
sampled_mig <- function(sampledloc,locations_sampled,gtreel,space, plot=0){
  ## create a matrix loc*loc with frequency of eachdim(space possible migration sampled
  migmat = matrix(0,nrow=dim(space)[2],ncol=dim(space)[2])
  ## get the distribution of locations for the root of the tree
  rootn = which(is.na(gtreel[,1]))
  rootlocs = rep(0,length(sampledloc[,1]))
  for (m in 1:length(sampledloc[,1])){ ## each sample step in MCMC
    history = treel2mig(gtreel,sampledloc[m,],space)
    for (n in 1:dim(history)[1]){
      migmat[history[n,1],history[n,2]] = migmat[history[n,1],history[n,2]] + 1
    }
    rootlocs[m] = sampledloc[m,rootn]
  }
  migmat = (migmat/length(sampledloc[,1]))*100
  ## build the matrix with the N most frequent sampled migrations (N number of simulated migrations)
  samp_migmat = matrix(0,nrow=dim(space)[2],ncol=dim(space)[2])
  num_sim_mig = (length(gtreel[,1])-1)/2 ## number of internal nodes
  samp_migmat[which(migmat>=sort(migmat,decreasing=TRUE)[num_sim_mig])] = 1
  ## build the matrix with the simulated migrations
  sim_migmat = matrix(0,nrow=dim(space)[2],ncol=dim(space)[2])
  history = treel2mig(gtreel,locations_sampled,space)
  for (n in 1:dim(history)[1]){
    sim_migmat[history[n,1],history[n,2]] = sim_migmat[history[n,1],history[n,2]] + 1
  }
  #corr_mig = cor(as.vector(samp_migmat),as.vector(sim_migmat))
  sim0sam0 = sum(sim_migmat==0 & samp_migmat==0) ## true negative
  sim0sam1 = sum(sim_migmat==0 & samp_migmat>0) ## false positive
  sim1sam1 = sum(sim_migmat>0 & samp_migmat>0) ## true positive
  sim1sam0 = sum(sim_migmat>0 & samp_migmat==0) ## false negative
  ## find the most frequent locations for the root
  surootlocs = sort(unique(rootlocs))
  rloccount = lapply(surootlocs,function(x) sum(rootlocs==x))
  rootnmod = surootlocs[sort(unlist(rloccount),decreasing=TRUE,index.return=TRUE)$ix[1]]
  ## average coordinates (green root location in simulation, red root location estimated)
  rootnave = c(sum(space[1,rootlocs])/length(rootlocs),sum(space[2,rootlocs])/length(rootlocs))
  if (plot==1){
    plot(space[1,],space[2,])
    points(space[1,locations_sampled[which(is.na(gtreel[,1])),1]],space[2,locations_sampled[which(is.na(gtreel[,1])),1]],col='green',pch=19)
    points(rootnave[1],rootnave[2],col='red',pch=24)
    points(space[1,rootnmod],space[2,rootnmod],col='red',pch='x')
  }
  return(list(sum(sim_migmat),sim1sam1,sim1sam0,sim0sam1,sim0sam0,round(space[,rootnmod],3),round(rootnave,3)))
}



## Do partial update in Tree by selecting first a depth and...
select_subtree <- function(treeStruct){
  runif(1,0,max(treeStruct$node_height))
  allow_update = (treeStruct$node_height<runif(1,0,max(treeStruct$node_height)) & treeStruct$node_height>0)
  return(allow_update)
}

##
### simulate history of migrations from the root with model allowing migration to occupied locations (proba lambda)
### gtreel: tree in format "louis"
### optional: a list of possible locations for each internal nodes and locations of tips
### getL: return likelihood of a specific set of locations
### model: each node is 2 dispersals potentially on father's loc (1) or exactly one dispersal out of father's loc (2) or one mig per node (4)
### samunif: if ==1, draw the locations of internal nodes uniformly in the list of possible locations, otherwise use the model
sim_history <- function(gtreel, space, sigma, lambda, tau, locations = 0, getL = 0, model = 1, samunif = 0, numtip = 50){
  
  ## to compute history likelihood
  history_proba_val = 0
  if (model==4 && locations[[1]][1]>0) { ## only one migration per internal node possible
    rootn = which(is.na(gtreel$nodes[,1]))
    locations_sampled = matrix(0,nrow=length(gtreel$nodes[,1]),ncol=(1+(dim(space)[1]))) ## location(destination) | distance(s)
    
    ## get the absolute height of each internal node
    inode_height = internal_node_height(gtreel)
    ## need to treat the nodes from the oldest to the most recent
    node_id = sort(inode_height,decreasing=TRUE,index.return=TRUE)$ix
    if (node_id[1]!=rootn){stop('sim_history(): error root')}
    
    ## keep track of colonised locations
    occupied = vector('numeric',ncol(space))
    if (length(locations[[rootn]])>1){ ## uniform to get the location at the root
      locations_sampled[rootn,1:2] = c( sample(locations[[rootn]],1), 0)
    }else{
      locations_sampled[rootn,1:2] = c( locations[[rootn]], 0)
    }
    occupied[locations_sampled[rootn,1]] = 1
    for (n in 1:length(node_id)) {
      if ( getL==1 || getL==2 ){ ploc = locations[[node_id[n]]]
      }else{ ploc = locations_sampled[node_id[n],1]} ## parent's location
      if ( getL==1 || getL==2 ){## return likelihood of a specific migration event out of the possible ones
        siblings = which(gtreel$nodes[,1]==node_id[n])
        if (length(siblings)==0) {next}
        if (getL==2 && length(which(gtreel$nodes[,1]==siblings[1]))==0 && length(which(gtreel$nodes[,1]==siblings[2]))==0) {next} ## hastings for loc sampling, ignore last migrations
        if (locations[[siblings[1]]]==ploc && locations[[siblings[2]]]!=ploc) {
          noddest = siblings[2]
          locdest = locations[[siblings[2]]]
        } else if (locations[[siblings[2]]]==ploc && locations[[siblings[1]]]!=ploc) {
          noddest = siblings[1]
          locdest = locations[[siblings[1]]]
        } else if (locations[[siblings[1]]]==ploc && locations[[siblings[2]]]==ploc) {
          noddest = siblings[1]
          locdest = locations[[siblings[1]]]
        } else if (locations[[siblings[1]]]!=ploc && locations[[siblings[2]]]!=ploc) { stop('sim_history(): error 3') }
        if (n==1){
          migration_event = new_migration( space, occupied, sigma, lambda, tau, 0, 0, c(ploc,locdest,0) )
        }else{
          migration_event = new_migration( space, occupied, sigma, lambda, tau, 0, 0, c(ploc,locdest,abs(inode_height[node_id[n]]-inode_height[node_id[n-1]])) )
        }
        history_proba_val = history_proba_val + log(migration_event[3])
        occupied[locdest] = occupied[locdest] + 1
        #print(migration_event)
        #print(history_proba_val)
        #print(occupied)
      }else{## sample a migration from the father's location if children nodes have more than 1 possible location, according to model (samunif!=1) or uniformly in the list of possible node location (samunif==1)
        migration_event = c(0,0,1,0,vector('numeric',dim(space)[1])) ## default: proba equals 1
        siblings = which(gtreel$nodes[,1]==node_id[n])
        loc_sib1 = locations[[siblings[1]]]
        loc_sib2 = locations[[siblings[2]]]
        if (length(loc_sib1)==1 && length(loc_sib2)==1){
          locations_sampled[siblings[1],1] = loc_sib1
          locations_sampled[siblings[2],1] = loc_sib2
          if (loc_sib1==ploc) {occupied[loc_sib2] = occupied[loc_sib2] + 1}
          else if (loc_sib2==ploc) {occupied[loc_sib1] = occupied[loc_sib1] + 1}
          else {stop('sim_history(): error 0')}
          ## which one is the parent's
        } else if (length(loc_sib1)>1 && length(loc_sib2)==1){
          if (loc_sib2==ploc) {
            if (samunif==1){ migration_event[2] = sample(loc_sib1,1) }
            else {
              migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib1, 0 )
            }
            locations_sampled[siblings[1],] = c(migration_event[2],migration_event[5:length(migration_event)])
            occupied[migration_event[2]] = occupied[migration_event[2]] +1
          }else{ ## loc_sib1 must contain ploc
            if (sum(loc_sib1==ploc)==0) stop('sim_history(): error 1')
            locations_sampled[siblings[1],1] = ploc
            occupied[loc_sib2] = occupied[loc_sib2] + 1
          }
          locations_sampled[siblings[2],1] = loc_sib2
        } else if (length(loc_sib1)==1 && length(loc_sib2)>1){
          if (loc_sib1==ploc) {
            if (samunif==1){ migration_event[2] = sample(loc_sib2,1) }
            else {
              migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib2)
            }
            locations_sampled[siblings[2],] = c(migration_event[2],migration_event[5:length(migration_event)])
            occupied[migration_event[2]] = occupied[migration_event[2]] + 1
          }else{ ## loc_sib2 must contain ploc
            if (sum(loc_sib2==ploc)==0) stop('sim_history(): error 2')
            locations_sampled[siblings[2],1] = ploc
            occupied[loc_sib1] = occupied[loc_sib1] + 1
          }
          locations_sampled[siblings[1],1] = loc_sib1
        } else if (sum(loc_sib1==ploc)>0 && sum(loc_sib2==ploc)>0){ ## draw which daughter stays according to number of leaves below that are father's location
          numer1 = sum(loc_sib1==ploc)
          numer2 = sum(loc_sib2==ploc)
          denom = numer1 + numer2
          sib_stays = sample(1:2,1,prob=c(numer1/denom,numer2/denom))
          if (sib_stays==1){
            locations_sampled[siblings[1],1] = ploc
            if (samunif==1){ migration_event[2] = sample(loc_sib2,1) }
            else {
              migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib2 , 0 )
           }
            locations_sampled[siblings[2],] = c(migration_event[2],migration_event[5:length(migration_event)])
            occupied[migration_event[2]] = occupied[migration_event[2]] + 1
          }else{
            locations_sampled[siblings[2],1] = ploc
            if (samunif==1){ migration_event[2] = sample(loc_sib1,1) }
            else {
              migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib1 , 0  )
              }
            locations_sampled[siblings[1],] = c(migration_event[2],migration_event[5:length(migration_event)])
            occupied[migration_event[2]] = occupied[migration_event[2]] + 1
          }
        } else if (sum(loc_sib1==ploc)==0 && sum(loc_sib2==ploc)>0){
          locations_sampled[siblings[2],1] = ploc
          if (samunif==1){ migration_event[2] = sample(loc_sib1,1) }
          else {
            migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib1, 0 )
          }
          locations_sampled[siblings[1],] = c(migration_event[2],migration_event[5:length(migration_event)])
          occupied[migration_event[2]] = occupied[migration_event[2]] + 1
        } else if (sum(loc_sib1==ploc)>0 && sum(loc_sib2==ploc)==0){
          locations_sampled[siblings[1],1] = ploc
          if (samunif==1){ migration_event[2] = sample(loc_sib2,1) }
          else {
            migration_event = new_migration( space, occupied, sigma, lambda, tau, ploc, loc_sib2, 0 )
          }
          locations_sampled[siblings[2],] = c(migration_event[2],migration_event[5:length(migration_event)])
          occupied[migration_event[2]] = occupied[migration_event[2]] + 1
        }
        history_proba_val = history_proba_val + log(migration_event[3])
      }
    }
    # calculate LogLikelihood
    history_proba_val = log(1/ncol(space)) + history_proba_val
    #cat(paste("\n",space[1,locations_sampled[rootn,1]],space[2,locations_sampled[rootn,1]],sigma[1],sigma[2],lambda,tau,history_proba_val,"\n"))
  } else { ## SIMULATION: simulate a sequence of migration events first, generate a tree second
    #rep_tlim = rep(0,100)
    #for (nrep in 1:100){
    #	rep_tlim[nrep] = sim_mig(space,sigma,0,tau,0)[[2]]
    #}
    #history = 0
    #while(is.vector(history)){ ## need at least 2 migrations, i.e. a matrix
    #	history = sim_mig(space,sigma,lambda,tau,quantile(rep_tlim,.95))[[1]]
    #}
    
    history = cbind(0)
    history = sim_mig(space,sigma,lambda,tau,0,numtip)[[1]] # limit to numtip leaves
    print(paste('number of simulated migration events:',nrow(history)-1))
    treeloc = mig2treel(history)
    gtreel = treeloc[[1]]
    #ici on a les nodes mais pas order
    locations_sampled = treeloc[[3]]
    # sort and reorder gtreel for consistent format
    reordered = sort_gtreel(gtreel, locations_sampled)
    gtreel = reordered[[1]]
    locations_sampled = gtreel$nodes.label
  }
  
  #print("gtreel")
  #print(gtreel)
  return(list(locations_sampled,history_proba_val,gtreel))
}

## simulate a sequence of migration events
sim_mig <- function(space, sigma, lambda, tau, timelimit, poplimit=0){
  occupied = vector('numeric', ncol(space))
  occupied[sample(1:ncol(space),1)] = 1 ## root ici pour choisir le racine

  history = matrix(0, 0, 4+dim(space)[1]) ## |dep|dest|proba|time2father|dist(|dist...)|
  time_elapsed = 0
  stop = 0
  while (stop==0){
    mig_event = new_migration(space, occupied, sigma, lambda, tau , 0 , 0 , 0)
    if(is.na(mig_event[4])){}
    time_elapsed = time_elapsed + mig_event[4]
    history = rbind(history, mig_event, deparse.level=0)
    occupied[mig_event[2]] = occupied[mig_event[2]] + 1
    if (timelimit > 0){
      if (time_elapsed > timelimit) {
        if (nrow(history) > 1){
          history = history[1:(nrow(history)-1),] ## remove last migration that goes above the limit
        }
        if (is.vector(history)){ ## control, at least 2 migrations needed (one is removed in mig2treel())
          return(list(history, time_elapsed))
        }
        history[nrow(history),4] = history[nrow(history),4] + (timelimit-sum(history[,4]))
        stop=1
      }
    } else if (poplimit>0) {
      if (sum(occupied)>=poplimit) {stop=1}
    } else {
      if (sum(occupied>0) == length(occupied)) {stop=1}
    }
  }
  if (nrow(history) == 1){
    history = as.vector(history)
  }
  return(list(history, time_elapsed))
}

##
## get some stats from MCMC output
stat_phyloland <- function(mcmco,plot_MCMC=0,save_MCMC=0){
  output_val = mcmco
  Nstep = output_val[length(output_val[,1]),1]
  if (plot_MCMC==1) {
    ## Attention aux variables non declarees
    #hist(output_val[,2],nclass=50,xlim=c(.1,3))
    hist(output_val[,2],nclass=50)
    par(mfrow=c(2,1))
    plot(1:length(output_val[,1]),output_val[,2],type="l")
    plot(output_val[,3],type="l",col="red")
    dev.copy2eps(file=paste("output_val",sigma_simul,"_",Nstep,".eps",sep=""))
  }
  if (save_MCMC==1) {
    filename = format(Sys.time(),"%a%d%b%Y_%H%M%S")
    ## need a first sample at 0 to be read correctly by tracer
    write.table(rbind(c(0,rep(0,Ndim),0,0,0),output_val),file=paste(filename,".log",sep=""), row.names=FALSE,col.names=c("Sample",paste("Sigma",1:Ndim),"Lambda","Tau","LogLikelihood"),sep="\t")
    print(paste(filename,".log"," written",sep=""))
    system(paste("tracer ",filename,".log ","&",sep=""))
  }
  stat = vector('numeric',dim(mcmco)[2]-4)
  smod = vector('numeric',dim(mcmco)[2]-4)
  slast = vector('numeric',dim(mcmco)[2]-4)
  smean = vector('numeric',dim(mcmco)[2]-4)
  squant = matrix(0,(dim(mcmco)[2]-4),7)
  ## SIGMAS ##
  for (n in 1:(dim(mcmco)[2]-4)){
    id = n+1
    ## p-value, number of samples below the threshold above which Ho can't be rejected
    #stat[n] = sum(mcmco[,id]>1.117) / nrow(mcmco)
    ## number of samples equal to 2.527
    stat[n] = sum(mcmco[,id]==2.527) / nrow(mcmco)
    ## ignore values equal to 2.527
    samplevals = mcmco[which(mcmco[,id]!=2.527),id]
    if (length(samplevals)>1){ ## no or only one value different to 2.527
      ## find the most frequent value (mode) of the distribution of sigma
      den = density(samplevals)
      smod[n] = den$x[which(den$y==max(den$y))[1]]
      ## quantiles
      squant[n,] = quantile(samplevals,probs=c(0,.05,.25,.5,.75,.95,1),names=FALSE)
      smean[n] = mean(samplevals)
      slast[n] = samplevals[length(samplevals)]
    }else if (length(samplevals)==1){ ## no or only one value different to 2.527
      ## find the most frequent value (mode) of the distribution of sigma
      smod[n] = samplevals
      ## quantiles
      squant[n,] = quantile(samplevals,probs=c(0,.05,.25,.5,.75,.95,1),names=FALSE)
      smean[n] = samplevals
      slast[n] = samplevals
    }else{
      smod[n] = 2.527
      squant[n,] = rep(2.527,7)
      smean[n] = 2.527
      slast[n] = 2.527
    }
    # 		## compute the area for the mode
    # 		norma_factor = pnorm(1,0,smodL)-pnorm(0,0,smodL)
    # 		x = sqrt( -2 * (smodL)^2 * log(norma_factor*sqrt(2*pi)*smodL) )
    # 		areaL = x + (pnorm(1,0,smodL)-pnorm(x,0,smodL))/norma_factor
  }
  ## LAMBDA ##
  ## find the most frequent value (mode) of the distribution of lambda
  den = density(mcmco[,(dim(mcmco)[2]-2)])
  lmod = den$x[which(den$y==max(den$y))[1]]
  lquant = quantile(mcmco[,(dim(mcmco)[2]-2)],probs=c(0,.05,.25,.5,.75,.95,1),names=FALSE)
  lmean = mean(mcmco[,(dim(mcmco)[2]-2)])
  llast = mcmco[nrow(mcmco),(dim(mcmco)[2]-2)]
  ## TAU ##
  ## find the most frequent value (mode) of the distribution of lambda
  den = density(mcmco[,(dim(mcmco)[2]-1)])
  tmod = den$x[which(den$y==max(den$y))[1]]
  tquant = quantile(mcmco[,(dim(mcmco)[2]-1)],probs=c(0,.05,.25,.5,.75,.95,1),names=FALSE)
  tmean = mean(mcmco[,(dim(mcmco)[2]-1)])
  tlast = mcmco[nrow(mcmco),(dim(mcmco)[2]-1)]
  
  ## FORMAT OUTPUT ##
  #strout = paste(pval[1],toString(squant[1,]),smod[1],smean[1],slast[1],sep=",")
  
  for (n in 1:(dim(mcmco)[2]-4)){
    if (n==1){ strout = paste(round(stat[n],3),toString(round(squant[n,],3)),round(smod[n],3),round(smean[n],3),round(slast[n],3),sep=",") }
    else if (n>1){ strout = paste(strout,round(stat[n],3),toString(round(squant[n,],3)),round(smod[n],3),round(smean[n],3),round(slast[n],3),sep=",") }
  }
  strout = paste(strout,toString(round(lquant,3)),round(lmod,3),round(lmean,3),round(llast,3),toString(round(tquant,3)),round(tmod,3),round(tmean,3),round(tlast,3),sep=",")
  return(strout)
}

### read an history of coalescent events (from most recent to oldest)
### write a tree file as Cayley format (parentheses) reading branch length
### ilabel==1: internal labels are coalescent tips, ==2: internal labels are in node_id
treestring_bl <- function(coalhistory,node_id=0,ilabel=1){ ## read node heights in coalhistory[,3]
#print(coalhistory)
#print(node_id)
	N = nrow(coalhistory)
	if (node_id[1]==0){node_id=rep("",1,N)}
	subtree = vector(mode="list", length=ceiling(N/2)) ## subtree strings
	subtreelength = vector(mode="numeric", length=ceiling(N/2)) ## subtree heights
	subtree[[1]] = paste("(t",coalhistory[1,1],":",coalhistory[1,3],",t",coalhistory[1,2],":",coalhistory[1,3],")",sep="")
	subtreelength[1] = coalhistory[1,3]
	numsubtree = 2 ## index for subtrees
	for (i in 2:N) {
#print("----------------------------")
#print(subtree)
#print(subtreelength)
		str = ""
		from = coalhistory[i,1]
		toward = coalhistory[i,2]
#print(paste(from,toward))
#print(subtreelength)
		## coalescent events before current node involving the one of the coalescing nodes
		fr = c(which(coalhistory[1:(i-1),1]==from), which(coalhistory[1:(i-1),2]==from))
		to = c(which(coalhistory[1:(i-1),1]==toward),which(coalhistory[1:(i-1),2]==toward))
		if (length(to)==0 & length(fr)==0) { ## two tips, create a new subtree string
			subtree[[numsubtree]] = paste("(t",coalhistory[i,1],":",coalhistory[i,3],",t",coalhistory[i,2],":",coalhistory[i,3],")",sep="")
			subtreelength[numsubtree] = coalhistory[i,3]
			numsubtree = numsubtree+1
		} else if (length(to)==0 & length(fr)>0) { ## toward is a tip
			## need to find the subtree containing "from"
			idsub = grep(paste("t",from,":",sep=""),subtree)
			if (ilabel==1){
				subtree[[idsub]] = paste("(",subtree[[idsub]],"t",from,":",coalhistory[i,3]-subtreelength[idsub],",t",toward,":",coalhistory[i,3],")",sep="")
			}else{
				id_node_idf = max(fr)
if (id_node_idf==58 || id_node_idf==68){
	browser()
}
				subtree[[idsub]] = paste("(",subtree[[idsub]],"t",node_id[id_node_idf],":",coalhistory[i,3]-subtreelength[idsub],",t",toward,":",coalhistory[i,3],")",sep="")
			}
			subtreelength[idsub] = coalhistory[i,3]
		} else if (length(to)>0 & length(fr)==0) { ## from is a tip
			## need to find the subtree containing "toward"
			idsub = grep(paste("t",toward,":",sep=""),subtree)
			if (ilabel==1){
				subtree[[idsub]] = paste("(",subtree[[idsub]],"t",toward,":",coalhistory[i,3]-subtreelength[idsub],",t",from,":",coalhistory[i,3],")",sep="")
			}else{
				id_node_idt = max(to)
if (id_node_idt==58 || id_node_idt==68){
	browser()
}
				subtree[[idsub]] = paste("(",subtree[[idsub]],"t",node_id[id_node_idt],":",coalhistory[i,3]-subtreelength[idsub],",t",from,":",coalhistory[i,3],")",sep="")
			}
			subtreelength[idsub] = coalhistory[i,3]
		} else { ## none is a tip
			idsubf = grep(paste("t",from,":",sep=""),subtree)
			idsubt = grep(paste("t",toward,":",sep=""),subtree)
			if (ilabel==1){
				subtree[[idsubf]] = paste("(",subtree[[idsubf]],"t",from,":",coalhistory[i,3]-subtreelength[idsubf],",",subtree[[idsubt]],"t",toward,":",coalhistory[i,3]-subtreelength[idsubt],")",sep="")
			}else{
				id_node_idf = max(fr)
				id_node_idt = max(to)
if (id_node_idf==58 || id_node_idf==68 || id_node_idt==58 || id_node_idt==68){
	browser()
}
				subtree[[idsubf]] = paste("(",subtree[[idsubf]],"t",node_id[id_node_idf],":",coalhistory[i,3]-subtreelength[idsubf],",",subtree[[idsubt]],"t",node_id[id_node_idt],":",coalhistory[i,3]-subtreelength[idsubt],")",sep="")
			}
			subtreelength[idsubf] = coalhistory[i,3]
			subtree[[idsubt]] = ""
			subtreelength[idsubt] = 0
		}
	}
#print(subtree)
	tree_par = subtree[[grep("t",subtree)]] ## only one subtree should remain, all the other should be empty
	if (ilabel==2){
		tree_par = paste(tree_par,"t",node_id[i],";",sep="")
	} else{
		tree_par = paste(tree_par,";",sep="")
	}
	return(tree_par)
}



##
### convert louis tree structure to the history of migrations
treel2mig <- function(gtreel, location, space){
  check_treel(gtreel, location)
  history = matrix(0, 0, 5+dim(space)[1]) ## |dep|dest|proba|time2father|dist(|dist...)|height
  rootn = which(is.na(gtreel$nodes[,1]))
  ## get the absolute height of each internal node
  inode_height = internal_node_height(gtreel)
  ## need to treat the nodes from the oldest to the most recent
  node_id = sort(inode_height, decreasing=TRUE, index.return=TRUE)$ix
  if (node_id[1] != rootn){stop('treel2mig error')}
  for (n in 1:length(node_id)) {
    ploc = location[[node_id[n]]]
    siblings = which(gtreel$nodes[,1] == node_id[n])
    if (length(siblings) == 0) {next}
    if (location[[siblings[1]]] == ploc && location[[siblings[2]]] != ploc) {
      noddest = siblings[2]
      locdest = location[[siblings[2]]]
    } else if (location[[siblings[2]]] == ploc && location[[siblings[1]]] != ploc) {
      noddest = siblings[1]
      locdest = location[[siblings[1]]]
    } else if (location[[siblings[1]]] == ploc && location[[siblings[2]]] == ploc) {
      noddest = siblings[1]
      locdest = location[[siblings[1]]]
    }
    distmig = as.vector(abs(space[,ploc]-space[,locdest]))
    history = rbind(history, c(ploc, locdest, NA, gtreel$nodes[noddest,4], distmig, inode_height[node_id[n]]))
  }
  return(history)
}

##
### convert louis tree structure to parenthesis (need library(ape))
treel2par <- function(gtreel){
	tree_phylo = lconverttree(gtreel)
	return(write.tree(tree_phylo))
}

##
### convert ape phylo tree structure (edges) to node index:parent|left child|right child|branch length to parent|
# 	gtreel = matrix(NA, (2*tree_phylo$Nnode + 1), 4)
# 	for (n in 1:(2*tree_phylo$Nnode)) {
# 		gtreel[tree_phylo$edge[n,2],1] = tree_phylo$edge[n,1]
# 		if (is.na(gtreel[tree_phylo$edge[n,1],2])) {
# 			gtreel[tree_phylo$edge[n,1],2] = tree_phylo$edge[n,2] }
# 		else {
# 			gtreel[tree_phylo$edge[n,1],3] = tree_phylo$edge[n,2] }
# 		gtreel[tree_phylo$edge[n,2],4] = tree_phylo$edge.length[n]
# 	}
# 	return(gtreel)



##
#numev=1000;ct=vector("numeric",numev);for(n in 1:numev){list_geneal=geneal(filecontent,numi);ct[n]=coal_time(list_geneal[[2]],1,2);}
#hist(ct[ct>0])




## test tree conversions
test_conv_tree <- function(plot=0){
	space_dim = 2
	space_size = 5
	sigma_simul = runif(space_dim,min=.1,max=2.527) ## see below for sigma upper boundaries
	lambda_simul = runif(1,min=.01,max=.99) ## .99
	tau_simul = runif(1,min=0.1,max=10)
	model = 4 ## dispersal model, full (1) or only one dispersal per node (2) or simulate migrations before tree (4)
	print(paste(' space_size=',space_size,' space_dim=',space_dim,' model=',model,
		' sigma=',toString(sigma_simul),' lambda=',lambda_simul,' tau=',tau_simul,sep=''))
	use_molseq = 0
	sim_data = sim_phyloland(0,0,0,0,space_size,sigma_simul,lambda_simul,tau_simul,use_molseq,space_dim,model)
	space = sim_data[[1]]
	gtreel = sim_data[[2]]
	locations_sim = sim_data[[3]]
	if (plot==1) {plotree(gtreel,locations_sim[,1])}
	# convert tree and back
	treepar = treel2par(gtreel)
	tree_phylo = read.tree(text=treepar)
	tipnames = as.numeric(gsub("[a-z]","",tree_phylo$tip))
	nodenames = as.numeric(gsub("[a-z]","",tree_phylo$node))
	if ( length(unique(tipnames))!=length(tipnames) ||
		length(unique(nodenames))!=length(nodenames) ){
		#nodenames[duplicated(nodenames)] #print duplicated references
		#print(gtreel)
		#print(tree_phylo$edge)
		return(list(gtreel,locations_sim,tree_phylo))
	}
	gtreel = converttree(tree_phylo)
	if (plot==1) {plotree(gtreel,locations_sim[,1])}
}

## show estimated locations from MCMC output
## return matrix with distance between true (sim) location vs estimated location for each node with height in the tree
showlocat <- function(sampledloc,gtreel,locations_sim,space,plotting=0){
	if (plotting==1){
		x11()
		plotree(gtreel,locations_sim[,1])
		x11()
	}
  # real_loc=????? ***BUG***
	gnode_height = internal_node_height(gtreel)
	if (dim(space)[1]==1) {## work in 1D
		sort_space_id = sort(space,decreasing=TRUE,index.return=TRUE)$ix
		xbar = barplot(hist(locations_sampled[,1],breaks=0:ncol(space),plot=FALSE)$counts[sort_space_id])
		title('estimated internal locations')
		x11();par(mfrow=c(2,1));
		for (node_id in (tnumi+1):dim(locations_sampled)[1]){
			tmp = histloc[(node_id-tnumi),]
			tmp = tmp[is.finite(tmp)]
			plot(1:length(tmp),space[tmp],type="l",ylim=c(0,1))
			points(vector('numeric',length(space)[1])*0,space) ## plot space
			points(0,space[real_loc[node_id]],col='red',pch=19) ## color real location
			#print(paste("int. node:",node_id,"est. acc. rate:",sum(diff(tmp)!=0)/(length(tmp)/(tnumi-1)))) ## on average, each internal nodes is treated Nstep/(tnumi-1) times
			#print(paste("int. node:",node_id,"est. acc. rate:",sum(diff(tmp)!=0)/sum(histdraw==node_id)))
			#print(paste("int. node:",node_id,"est. acc. rate:",sum(acc_loc[node_id-tnumi,])/100))
			#hist(tmp,breaks=0:ncol(space))
			sort_space_id = sort(space,decreasing=TRUE,index.return=TRUE)$ix
			xbar = barplot(hist(tmp,breaks=0:ncol(space),plot=FALSE)$counts[sort_space_id])
			points(xbar[which(sort_space_id==real_loc[node_id])],0,col='red',pch=19) ## color real location
			title(paste("int. node:",node_id,"depth",gnode_height[node_id]/max(gnode_height)))
			#Sys.sleep(.5)
			readline()
		}
	}else if (dim(space)[1]==2){
		height_dist = matrix(0,length(sampledloc[1,]),2)
		for (n in 1:length(sampledloc[1,])){ ## each node of the tree
			if (plotting==1){
				plot(space[1,],space[2,],pch='.')
			}
			locn = vector('numeric',length(sampledloc[,1]))
			for (m in 1:length(sampledloc[,1])){ ## each sample of node locations
				locn[m] = sampledloc[m,n]
			}
			uniklocn = unique(locn)
			occulocn = vector('numeric',length(uniklocn))
			for (l in 1:length(uniklocn)){
				occulocn[l] = sum(locn==uniklocn[l])
				print("occ")
				print(occulocn)
				print(occulocn[l])
				print("cex")
				print((1+(10*occulocn[l])/length(sampledloc[,1])))
				print("len")
				print(length(sampledloc[,1]))
				if (plotting==1){
					points(space[1,uniklocn[l]],space[2,uniklocn[l]],cex=(1+(10*occulocn[l])/length(sampledloc[,1])))
					text(space[1,uniklocn[l]],space[2,uniklocn[l]],uniklocn[l])
				}
			}
			if (plotting==1){
				points(space[1,locations_sim[n,1]],space[2,locations_sim[n,1]],col='green',pch=19)
			}
			#title(paste('node ',n,', height ',gnode_height[n],' (',gnode_height[which(is.na(gtreel[,1]))],')',sep=''))
			height_dist[n,1] = round(100*gnode_height[n]/gnode_height[which(is.na(gtreel[,1]))],3)
			#estim = uniklocn[which(occulocn==max(occulocn))[1]]
			estim = uniklocn[match(max(occulocn),occulocn)] # in case of ties, take the first one
			print(estim)
			height_dist[n,2] = sqrt( (space[1,locations_sim[n,1]]-space[1,estim])^2 + 
						(space[2,locations_sim[n,1]]-space[2,estim])^2 )
			if (plotting==1){
				title(paste('node ',n,', height ',height_dist[n,1],'% distance ',round(height_dist[n,2],3),sep=''))
				readline()
			}
		}
	}
	return(height_dist)
}

## simulate landscape and return random pairwise Euclidean distances
#### ? data #####
sim_landscape_dist <- function(testcorr=0,printfile=0){
	pwdist = vector('numeric',1000)
	for (i in 1:1000){
		space_size = round(runif(1,min=5,max=100))
		space = create_landscape(space_size,2)
		loc1 = sample(space_size,1)
		loc2 = sample(space_size,1)
		pwdist[i] = sqrt((space[1,loc1]-space[1,loc2])^2+(space[2,loc1]-space[2,loc2])^2)
	}
	dmod = sqrt((data[,61]-data[,63])^2+(data[,62]-data[,64])^2)
	d1 = density(pwdist)
	d2 = density(dmod)
	plot(d1,type='l',main="",xlab="",ylab="",lwd=2,yaxt="n",xlim=c(-.2,1.5),col='grey80')
	polygon( c(d1$x[1],d1$x,d1$x[length(d1$x)]), c(0,d1$y,0), col="grey80", lty=0) #color the area below the curve
	points(d2,type='l',lwd=2)
	if (printfile==1){
		dev.print(file="/home/louis/Documents/FRDF_steph_2010/Simulations/Retrieved/9Feb2012/root_distance.eps",height=10)
	}
	if (testcorr==1){
		n1 = 1000
		n2 = 1000
		sprintf("obs. correlation root mean coord1: %.2f",cor(data[,61],data[,65]))
		sprintf("obs. correlation root mean coord2: %.2f",cor(data[,62],data[,66]))
		corrdist = vector('numeric',n1)
		for (j in 1:n1){
			coordrnd = matrix('numeric',n2,2)
			for (i in 1:n2){
				space_size = round(runif(1,min=5,max=100))
				space = create_landscape(space_size,2)
				loc1 = sample(space_size,1)
				loc2 = sample(space_size,1)
				coordrnd[i,] = c(space[1,loc1],space[2,loc2])
			}
			corrdist[j] = cor(coordrnd[,1],coordrnd[,2])
		}
		d3 = density(corrdist)
		plot(d3,type='l',main="",xlab="",ylab="",lwd=2,yaxt="n",xlim=c(-.6,.6),col='grey80')
		polygon( c(d3$x[1],d3$x,d3$x[length(d3$x)]), c(0,d3$y,0), col="grey80", lty=0) #color the area below the curve
		abline(v=cor(data[,61],data[,65]))
		abline(v=cor(data[,62],data[,66]))
		if (printfile==2){
			dev.print(file="/home/louis/Documents/FRDF_steph_2010/Simulations/Retrieved/9Feb2012/root_distance_corr.eps",height=10)
		}
	}
}





	

	
	

##### Return new edges length gtreel->tree #####
scale_tree <- function(tree_phylo,M=0.5,alpha=1,save_data=1,file_save=NA){
	lengths=tree_phylo$edge.length 
	lengths_mem=lengths                                   
	T<-height(which(is.na(converttree(tree_phylo)[,1])),converttree(tree_phylo)) ## hauteur de la racine
	#M=1 ## nombre de mutations entre une feuille et la racine
	#a=runif(1,0.9,4)
	c=vector('numeric',0)
	for (i in 1:length(lengths)){
		c[i]=rgamma(1,shape=alpha)
		lengths[i]=((lengths[i]*M*c[i])/T)
		tree_phylo$edge.length[i]=lengths[i]
		}
	if (save_data==1){
	liste=list(lengths_mem,c)
	write.table(liste, file=file_save,append = FALSE, quote = FALSE, sep = " ", na = "-", col.names=c("length","c"), row.names = FALSE,eol='\n')
		}
	return(list(tree_phylo,c))
	}

	
##### Return new edges length tree->gtreel #####
unscale_tree <- function(tree_phylo,gtreel){
	height_root_treel<-height(which(is.na(gtreel[,1])),gtreel)
	height_root_phylo<-height(which(is.na(converttree(tree_phylo)[,1])),converttree(tree_phylo))
	new_tree_phylo<-tree_phylo
	new_tree_phylo$edge.length=((tree_phylo$edge.length*height_root_treel)/(height_root_phylo))
	return(new_tree_phylo)
	}



##### Save tree #####
save_tree <- function(gtreel,fileOUT){
	write.table(gtreel, file = fileOUT, append = FALSE, quote = FALSE, sep = " ", na = "-", col.names=FALSE, row.names = FALSE,eol='')
	}

##### Save gtreel #####
save_treel <- function(gtreel,file1){
	write.table(gtreel,file = file1, append = FALSE, quote = FALSE, sep = " ", na = "-", col.names=FALSE, row.names = FALSE,eol='\n')
	}
	
##### Save locations #####
save_loc <- function(locations,file1){
	write.table(locations,file = file1, append = FALSE, quote = FALSE, sep = " ", na = "-", col.names=FALSE, row.names = FALSE,eol='\n')
	}




	
##### Return commande line for SeqGen #####
seq_gen <- function(m='HKY',z=600,z_alea=1,l=1000,t=4,f=rep(0.25,4),nexus=T,freq_alea=1,fileIN,fileOUT){
	if (z_alea==1){					
		z=sample(1:100000000,1)
		}
	if (freq_alea==1){
		f=freq_nucl(save_freq=0)
		}
	if (nexus){
		line=paste('./seq-gen -m',m,' -z',z,' -l',l,' -t',t,' -f ',f[1],' ',f[2],' ',f[3],' ',f[4],' -on',' < ',fileIN,' > ',fileOUT,sep="",eol='')
		}
	else{
		line=paste('./seq-gen -m',m,' -z',z,' -l',l,' -t',t,' -f ',f[1],' ',f[2],' ',f[3],' ',f[4],' < ',fileIN,' > ',fileOUT,sep="",eol='')
		}
	return(line)
	}
	
###### Return commande line for SeqGen : grid version #####
#seq_gen <- function(m='HKY',z=600,z_alea=1,l=1000,t=4,f=rep(0.25,4),nexus=T,freq_alea=1,fileIN,fileOUT){
#	if (z_alea==1){					
#		z=sample(1:100000000,1)
#		}
#	if (freq_alea==1){
#		f=freq_nucl(save_freq=0)
#		}
#	if (nexus){
#		line=paste('/share/apps/seqgen/seq-gen -m',m,' -z',z,' -l',l,' -t',t,' -f ',f[1],' ',f[2],' ',f[3],' ',f[4],' -on',' < ',fileIN,' > ',fileOUT,sep="",eol='')
#		}
#	else{
#		line=paste('/share/apps/seqgen/seq-gen -m',m,' -z',z,' -l',l,' -t',t,' -f ',f[1],' ',f[2],' ',f[3],' ',f[4],' < ',fileIN,' > ',fileOUT,sep="",eol='')
#		}
#	return(line)
#	}
	
	
##### Return file xml for BEAST #####
to_xml <- function(seqbase_xml,seqbase2_xml,tree_par,seq_nex,fileXML,topo_given=F,id_filena){
	
	if(topo_given==T){
		Seq_Base=readLines(seqbase_xml)
		seq=readLines(seq_nex)
		new_xml=Seq_Base
	
		a=(which(seq=="\tMatrix")+1)
    	z=(which(seq=="\t;")-1)
    	nb_tips=z-a+1
    
    	tax=vector('character',nb_tips)
    	sequences=vector('character',nb_tips)
    	for (i in 1:nb_tips){
			s=strsplit(seq[i+a-1],split=" ")[[1]]
			tax[i]=s[1]
    		if (s[2]==""){
    			s=strsplit(seq[i+a-1],split=" ")[[1]][-which(strsplit(seq[i+a-1],split=" ")[[1]]=="")]
				}
			sequences[i]=s[2]
			}

    	ind_nb_tax1=which(Seq_Base=="\t<!-- ntax=12                                                                 -->")
    	start_tax=c((which(Seq_Base=="\t<taxa id=\"taxa\">"))+1,which(Seq_Base=="\t<taxa id=\"untitled1\">")+1)
    	stop_tax=which(Seq_Base=="\t</taxa>")-1
    	nb_tax=stop_tax[1]-start_tax[1]+1
    	ind_nb_tax2=which(Seq_Base=="\t<!-- ntax=12 nchar=1000                                                      -->")   
    	start_seq=which(Seq_Base=="\t<alignment id=\"alignment\" dataType=\"nucleotide\">")+1
    	stop_seq=which(Seq_Base=="\t</alignment>")-1
    	nb_seq=stop_seq[1]-start_seq[1]+1 
    	start_tree=which(Seq_Base== "\t<!-- The user-specified starting tree in a newick tree format.               -->")
    	stop_tree=which(Seq_Base=="\t</newick>")
    
    	new_xml=new_xml[c(-(start_tree:stop_tree),-(start_seq:stop_seq),-(ind_nb_tax2),-(start_tax[2]:stop_tax[2]),-(start_tax[1]:stop_tax[1]),-(ind_nb_tax1))]

		insert_tree=vector("character",4)
    	insert_tree[1] = "\t<!-- The user-specified starting tree in a newick tree format.               -->"
    	insert_tree[2] = "\t<newick id=\"startingTree\" usingDates=\"false\">"  
    	insert_tree[3] = paste("\t\t",tree_par,sep="")
    	insert_tree[4] = "\t</newick>"
    	
    	nb_taxa1=paste("\t<!-- ntax=",nb_tips,"                                                                 -->",sep="")
    	taxa1=vector('character',nb_tips)
    	taxa2=vector('character',nb_tips)
    	for (i in 1:nb_tips){
    		taxa1[i]=paste("\t\t<taxon id=\"",tax[i],"\"/>",sep="")
    		taxa2[i]=paste("\t\t<taxon idref=\"",sort(tax)[i],"\"/>",sep="")
    		}
    	nb_taxa2=paste("\t<!-- ntax=",nb_tips," nchar=1000                                                      -->",sep="")
    	insert_seq=vector('character',(4*nb_tips))
    	for (i in 0:(nb_tips-1)){
    		insert_seq[4*i+1]="\t\t<sequence>"   
        	insert_seq[4*i+2]=paste("\t\t\t<taxon idref=\"",tax[i+1],"\"/>",sep="")
        	insert_seq[4*i+3]=paste("\t\t\t",sequences[i+1],sep="")
        	insert_seq[4*i+4]="\t\t</sequence>"
    		}
			
		start_tax=c((which(new_xml=="\t<taxa id=\"taxa\">")),which(new_xml=="\t<taxa id=\"untitled1\">"))
    	start_seq=which(new_xml=="\t<alignment id=\"alignment\" dataType=\"nucleotide\">")
    	ind_tree=which(new_xml=="\t</constantSize>")+1

    	ind_nb_tax1=start_tax[1]-1
    	ind_nb_tax2=start_seq-1
	
		new_xml=insert(x=(insert(x=(insert(x=(insert(x=(insert(x=insert(x=new_xml,ats=(ind_tree+1),values=insert_tree),ats=(start_seq+1),values=insert_seq)),ats=(ind_nb_tax2+1),values=nb_taxa2)),ats=(start_tax[2]+1),values=taxa2)),ats=(start_tax[1]+1),values=taxa1)),ats=(ind_nb_tax1+1),values=nb_taxa1)
	
	
	
		ind_line_names=vector("numeric",3)
		ind_line_names[1]=which(new_xml=="\t<mcmc id=\"mcmc\" chainLength=\"10000000\" autoOptimize=\"true\" operatorAnalysis=\"untitled.ops\">" )
		ind_line_names[2]=which(new_xml=="\t\t<log id=\"fileLog\" logEvery=\"1000\" fileName=\"untitled.log\" overwrite=\"false\">")
		ind_line_names[3]=which(new_xml=="\t\t<logTree id=\"treeFileLog\" logEvery=\"1000\" nexusFormat=\"true\" fileName=\"untitled.trees\" sortTranslationTable=\"true\">" )	
		line_names=vector("character",3)
		line_names[1]=paste("\t<mcmc id=\"mcmc\" chainLength=\"10000000\" autoOptimize=\"true\" operatorAnalysis=\"","trees1_",id_filena,".ops","\">",sep="") 
		line_names[2]=paste("\t\t<log id=\"fileLog\" logEvery=\"1000\" fileName=\"","trees1_",id_filena,".log","\" overwrite=\"false\">",sep="")
		line_names[3]=paste("\t\t<logTree id=\"treeFileLog\" logEvery=\"1000\" nexusFormat=\"true\" fileName=\"","trees1_",id_filena,".trees","\" sortTranslationTable=\"true\">",sep="")
		
		
		new_xml=new_xml[c(-(ind_line_names[1]),-(ind_line_names[2]),-(ind_line_names[3]))]
		new_xml=insert(x=(insert(x=insert(x=new_xml,ats=ind_line_names[1],values=line_names[1]),ats=ind_line_names[2],values=line_names[2])),ats=ind_line_names[3],values=line_names[3])
	
		writeLines(new_xml,fileXML)
		}
		
		if (topo_given==F){
		Seq_Base=readLines(seqbase2_xml)
		seq=readLines(seq_nex)
		new_xml=Seq_Base
	
		a=(which(seq=="\tMatrix")+1)
    	z=(which(seq=="\t;")-1)
    	nb_tips=z-a+1
    
    	tax=vector('character',nb_tips)
    	sequences=vector('character',nb_tips)
    	for (i in 1:nb_tips){
			s=strsplit(seq[i+a-1],split=" ")[[1]]
			tax[i]=s[1]
    		if (s[2]==""){
    			s=strsplit(seq[i+a-1],split=" ")[[1]][-which(strsplit(seq[i+a-1],split=" ")[[1]]=="")]
				}
			sequences[i]=s[2]
			}

    	ind_nb_tax1=which(Seq_Base=="\t<!-- ntax=12                                                                 -->")
    	start_tax=c((which(Seq_Base=="\t<taxa id=\"taxa\">"))+1,which(Seq_Base=="\t<taxa id=\"untitled1\">")+1)
    	stop_tax=which(Seq_Base=="\t</taxa>")-1
    	nb_tax=stop_tax[1]-start_tax[1]+1
    	ind_nb_tax2=which(Seq_Base=="\t<!-- ntax=12 nchar=1000                                                      -->")   
    	start_seq=which(Seq_Base=="\t<alignment id=\"alignment\" dataType=\"nucleotide\">")+1
    	stop_seq=which(Seq_Base=="\t</alignment>")-1
    	nb_seq=stop_seq[1]-start_seq[1]+1 
    
    	new_xml=new_xml[c(-(start_seq:stop_seq),-(ind_nb_tax2),-(start_tax[2]:stop_tax[2]),-(start_tax[1]:stop_tax[1]),-(ind_nb_tax1))]
  	
    	nb_taxa1=paste("\t<!-- ntax=",nb_tips,"                                                                 -->",sep="")
    	taxa1=vector('character',nb_tips)
    	taxa2=vector('character',nb_tips)
    	for (i in 1:nb_tips){
    		taxa1[i]=paste("\t\t<taxon id=\"",tax[i],"\"/>",sep="")
    		taxa2[i]=paste("\t\t<taxon idref=\"",sort(tax)[i],"\"/>",sep="")
    		}
    	nb_taxa2=paste("\t<!-- ntax=",nb_tips," nchar=1000                                                      -->",sep="")
    	insert_seq=vector('character',(4*nb_tips))
    	for (i in 0:(nb_tips-1)){
    		insert_seq[4*i+1]="\t\t<sequence>"   
        	insert_seq[4*i+2]=paste("\t\t\t<taxon idref=\"",tax[i+1],"\"/>",sep="")
        	insert_seq[4*i+3]=paste("\t\t\t",sequences[i+1],sep="")
        	insert_seq[4*i+4]="\t\t</sequence>"
    		}
			
		start_tax=c((which(new_xml=="\t<taxa id=\"taxa\">")),which(new_xml=="\t<taxa id=\"untitled1\">"))
    	start_seq=which(new_xml=="\t<alignment id=\"alignment\" dataType=\"nucleotide\">")
    	ind_nb_tax1=start_tax[1]-1
    	ind_nb_tax2=start_seq-1
	
		new_xml=insert(x=(insert(x=(insert(x=(insert(x=(insert(x=new_xml,ats=(start_seq+1),values=insert_seq)),ats=(ind_nb_tax2+1),values=nb_taxa2)),ats=(start_tax[2]+1),values=taxa2)),ats=(start_tax[1]+1),values=taxa1)),ats=(ind_nb_tax1+1),values=nb_taxa1)
	
	
		ind_line_names=vector("numeric",3)
		ind_line_names[1]=which(new_xml=="\t<mcmc id=\"mcmc\" chainLength=\"10000000\" autoOptimize=\"true\" operatorAnalysis=\"untitled.ops\">" )
		ind_line_names[2]=which(new_xml=="\t\t<log id=\"fileLog\" logEvery=\"1000\" fileName=\"untitled.log\" overwrite=\"false\">")
		ind_line_names[3]=which(new_xml=="\t\t<logTree id=\"treeFileLog\" logEvery=\"1000\" nexusFormat=\"true\" fileName=\"untitled.trees\" sortTranslationTable=\"true\">" )	
		line_names=vector("character",3)
		line_names[1]=paste("\t<mcmc id=\"mcmc\" chainLength=\"10000000\" autoOptimize=\"true\" operatorAnalysis=\"","trees2_",id_filena,".ops","\">",sep="") 
		line_names[2]=paste("\t\t<log id=\"fileLog\" logEvery=\"1000\" fileName=\"","trees2_",id_filena,".log","\" overwrite=\"false\">",sep="")
		line_names[3]=paste("\t\t<logTree id=\"treeFileLog\" logEvery=\"1000\" nexusFormat=\"true\" fileName=\"","trees2_",id_filena,".trees","\" sortTranslationTable=\"true\">",sep="")
		
		
		new_xml=new_xml[c(-(ind_line_names[1]),-(ind_line_names[2]),-(ind_line_names[3]))]
		new_xml=insert(x=(insert(x=insert(x=new_xml,ats=ind_line_names[1],values=line_names[1]),ats=ind_line_names[2],values=line_names[2])),ats=ind_line_names[3],values=line_names[3])
	
		writeLines(new_xml,fileXML)
		}

	
	}
	






# Checks if the tree is contained in the list of trees already sampled
test_tree <- function(tree_phylo, list_trees, list_nb, list_ind, n, freq){
	equal <- lapply(list_trees, "all.equal.phylo", current = tree_phylo, use.edge.length = FALSE)
	if (length(which(equal==TRUE)) > 0){
		list_nb[which(equal == TRUE)] = list_nb[which(equal == TRUE)] + 1
		list_ind[[which(equal == TRUE)]][length(list_ind[[which(equal == TRUE)]]) + 1] = n/freq
	}else{
		list_trees[[length(list_trees) + 1]] <- tree_phylo
		list_nb[length(list_nb) + 1] = 1	
		for ( i in 1:length(list_ind) ){
			if (length(list_ind[[i]]) == 0){
				list_ind[[i]][1] = n/freq
				break	
			}	
		}
	}
	return(list(list_trees, list_nb, list_ind))
}
	
# likelihood extraction
treeLikeli <- function(line, pattern){
	if(regexpr(pattern,line)[1] < 0 ){
		stop(" in interface : wrong 'pattern_trees_likelihood'",call.=FALSE)
	}
	start = regexpr(pattern,line)[1] + nchar(pattern)
	stop = regexpr("]",line)[1] -1
	Lchar = strsplit(line, split = "")[[1]][start:stop]
	return(as.numeric(paste(Lchar, collapse="")))
}

# compute pairwise distances of space
space_dist <- function(space, dmethod="euclidean"){
  mat_Dists = list()
  if (dmethod=="distkm"){ # separate latitudinal and longitudinal kernels by average other dimension when computing earth distance in km
    mat_Dists[[1]]=array(0,c(ncol(space),ncol(space)))
    mat_Dists[[2]]=array(0,c(ncol(space),ncol(space)))
    for (a in 1:(ncol(space)-1)){
      for (b in (a+1):ncol(space)){
        mat_Dists[[1]][a,b] = distkm( space[1,a], space[1,b], (space[2,a]+space[2,b])/2, (space[2,a]+space[2,b])/2 ) 
        mat_Dists[[1]][b,a] = mat_Dists[[1]][a,b]
        # pa=space[,a]
        # pb=space[,b]
        # pa=rev(pa)
        # pb=rev(pb)
        # angle=bearing(pa,pb)
        # mi_dist=distGeo(pa,pb)/2
        # PC=destPoint(pa,angle,mi_dist)
        # longc=PC[1]
        # latc=PC[2]
        # longa=space[2,a]
        # longb=space[2,b]
        # p1=c(longb,latc)
        # p2=c(longa,latc)
        # dist_mean=distGeo(p1,p2)/1000
        # mat_Dists[[2]][a,b]=dist_mean
        mat_Dists[[2]][a,b] = distkm( (space[1,a]+space[1,b])/2, (space[1,a]+space[1,b])/2, space[2,a], space[2,b] )
        mat_Dists[[2]][b,a] = mat_Dists[[2]][a,b]
      }
    }
    
  }else{
    for (n in 1:dim(space)[1]){
      mat_Dists[[n]] = dist(space[n,], diag=T, upper=T, method=dmethod)
    }
  }
  return(mat_Dists)
}

# find upper limit for sigma so that area of normal density represent 99% of the uniform U(0,dmax)
sigma_upperlim <- function(dmax,aire=.95,n_itera=100,sigma_min=1e-3,sigma_max=1e3,tries=1e4,plotting=0){
  for (i in 1:n_itera){
    testsigma = seq(sigma_min,sigma_max,length.out=tries)
    proba1 = matrix(0,nrow=length(testsigma),ncol=2)
    #      proba1[i,2]=dnorm(dmin,0,testsigma[i])/dnorm(0,0,testsigma[i])
    for(i in 1:length(testsigma)){
      proba1[i,1] = testsigma[i]
      #norma_factor = pnorm(0,0,sqrt(testsigma[i]))
      norma_factor = pnorm(dmax,0,sqrt(testsigma[i]))-pnorm(0,0,sqrt(testsigma[i]))
      #if ( (norma_factor * sqrt(2*pi) * sqrt(testsigma[i])) > 1 ) {print(paste(testsigma[i],norma_factor,"\n"));break;}
      x = sqrt( -2 * testsigma[i] * log(norma_factor * (1/dmax) * sqrt(2*pi) * sqrt(testsigma[i])) ) # +0 (because the mean equals to 0!)
      proba1[i,2] = x*(1/dmax) + (pnorm(dmax,0,sqrt(testsigma[i]))-pnorm(x,0,sqrt(testsigma[i])))/norma_factor
    }
    # error control if need to scale
    #return(proba1)
    #print(paste("proba1",proba1[length(testsigma),2]))
    if (proba1[length(testsigma),2]<aire) {
      #stop("sigma_upperlim(): distances are too large, need to be scaled")
      sigma_min = sigma_max
      sigma_max = sigma_max^2
    }else if (proba1[1,2]>aire) {
      #stop("sigma_upperlim(): distances are too short, need to be scaled")
      sigma_max = sigma_min
      sigma_min = sigma_min/2
    }else{
      # upper limit
      n = which(proba1[,2]>=aire)[1]
      sigma_min = proba1[n-1,1]
      sigma_max = proba1[n,1]
    }
  }
  # upper limit
  sig_upper_limit = testsigma[which(proba1[,2]>=aire)][1]
  if (plotting==1){
    xdim=seq(0,dmax*1.1,length.out=100)
    plot(xdim,dnorm(xdim,0,.3)/(pnorm(1,0,.3)-pnorm(0,0,.3)),type='l',lwd=2,xlab="distance",ylab="probability density function")
    abline(v=dmax)
    abline(h=1)
    points(xdim,dnorm(xdim,0,sig_upper_limit)/(pnorm(1,0,sig_upper_limit)-pnorm(0,0,sig_upper_limit)),type="l",lwd=2,col="red")
    text(c(dmax,dmax), c(dnorm(dmax*1.1,0,.3)/(pnorm(1,0,.3)-pnorm(0,0,.3)),dnorm(dmax*1.1,0,sig_upper_limit)/(pnorm(1,0,sig_upper_limit)-pnorm(0,0,sig_upper_limit))),
         labels = c("0.3",sig_upper_limit))
  }
  return(sig_upper_limit)
  #return(proba1)
}



# write simulation data to text files
# simulated_history is the output of sim_history()
write_sim_history <- function(simulated_history,space,file_tree="tree.tre",file_location="locations.txt"){
  library("stringr")
  library("gsubfn")
  library("ape")
  locations_sim = simulated_history[[1]]
  gtreel = simulated_history[[3]]
  tree_par = treel2par(gtreel)
  tree_phylo = read.tree(text=tree_par)
  tips = tree_phylo$tip
  tips = gsub("t","",tips)
  tips_location = locations_sim[,1][tips] # replace node number by location index
  #node_list = as.numeric(unlist(lapply(str_extract_all(tree_par,"t[0-9]*[:|;]"),function(x) gsub("[t:;]","",x))))
  #node_location = locations_sim[,1][node_list]
  tree_par = gsubfn('t([0-9]*):',function(x) paste(as.numeric(x),'-L',locations_sim[,1][as.numeric(x)],':',sep=''),tree_par)
  tree_par = gsubfn('t([0-9]*);',function(x) paste(as.numeric(x),'-L',locations_sim[,1][as.numeric(x)],';',sep=''),tree_par)
  cat(tree_par,file=file_tree)
  #renamed_node_list = unlist(lapply(str_extract_all(tree_par,"[)(,][0-9]*-L[0-9]*[:|;]"),function(x) gsub("[(),:;]","",x)))
  #write.table( cbind(renamed_node_list,t(space[,node_location])), file=file_location, sep=" ", col.names=F, row.names=F, quote=F )
  # only write tips locations
  tree_phylo2 = read.tree(text=tree_par)
  renamed_tips = tree_phylo$tip
  # NEED TO CONTROL THAT tips AND renamed_tips MAP I.E. THAT TIPS ARE ALWAYS LISTED IN THE SAME ORDER BY read.tree...
  write.table( cbind(renamed_tips,t(space[,tips_location])), file=file_location, sep=" ", col.names=F, row.names=F, quote=F )
}

# caculate tree height at which the space is fully colonised for the first time
saturation_height <- function(history,space){
  occupied = rep(0, ncol(space))
  height_id = length(history[1,]) # variable with number of dimensions
  for (n in 1:length(history[,height_id])){
    if (n==1){ # root
      occupied[history[n,1]] = occupied[history[n,1]] + 1
    }
    occupied[history[n,2]] = occupied[history[n,2]] + 1
    if (sum(occupied==0)==0){
      return(history[n,height_id])
    }
  }
  return(0)
}


# sample a set of tips in gtreel and locations_sampled
sample_tips <- function(gtreel, locations_sampled, space, ntip){
  tips = which(is.na(gtreel[,2]))
  if (ntip>length(tips)) stop('sample_tips: number of tips to sample is greater than total number of tips in the initial tree')
  history = treel2mig(gtreel, locations_sampled, space)
  while ((nrow(history)+1)>ntip){
    inode_height = internal_node_height(gtreel)
    node_height_s = sort(inode_height,index.return=T,decreasing=T)
    history = treel2mig(gtreel, locations_sampled, space)
    # select a tip at random
    tip_del = sample(which(is.na(gtreel[,2])),1)
    # get the rank in history (node height) of the parent that gave birth to the tip
    mig_id_del = which(node_height_s$ix==gtreel[tip_del,1])
    # discard the corresponding migration event and all subsequent migrations starting from it
    #print(paste(">",mig_id_del))
    history = delete_mig_history(history,mig_id_del)
    treeloc = mig2treel(history)
    gtreel = treeloc[[1]]
    locations_sampled = treeloc[[3]]
  }
  ordered = sort_gtreel(gtreel,locations_sampled)
  gtreel = ordered[[1]]
  locations_sampled = ordered[[2]]
  check_treel(gtreel,locations_sampled)
  return(list(gtreel, locations_sampled))
}



# sort gtreel so that tips come first, then root, then internal nodes
sort_gtreel <- function(gtreel, locations){
  gtreelO = list(nodes = NA, nodes.label=NA)
  gtreelO$nodes = matrix(0,nrow=nrow(gtreel$nodes),ncol=4)
  ntips = ( length(gtreel$nodes[,1])+1 )/2
  nnodes = ntips-1
  # save original indexes
  idtips = which(is.na(gtreel$nodes[,2]))
  rootn = which(is.na(gtreel$nodes[,1]))
  idintnodes = which(is.na(gtreel$nodes[,2])==F & is.na(gtreel$nodes[,1])==F) # all internal nodes
  gtreelO$nodes[1:ntips,2:3] = cbind(rep(NA,ntips),rep(NA,ntips)) # tips first (no children)
  gtreelO$nodes[ntips+1,c(1,4)] = NA # then root (no parent)
  gtreelO$nodes[1:ntips,4] = gtreel$nodes[idtips,4] # height for tips stay the same
  
  ## get the absolute height of each internal node
  inode_height = internal_node_height(gtreel)
  intnode_ordered = idintnodes[sort(inode_height[idintnodes],index.return=T)$ix] # internal nodes from most recent to oldest
  mapintnode = cbind( c(rootn,intnode_ordered), (ntips+1):(ntips+nnodes)) # keep the mapping between old and new indexes, including the root
  idintnodes0 = ntips+1+1
  for (n in intnode_ordered) { # go through internal nodes, from most recent to oldest
    if(is.na(gtreel$nodes[gtreel$nodes[n,2],2])){ # internal nodes of a tip in the original tree
      filsg = which(idtips==gtreel$nodes[n,2])
    }else{
      filsg = mapintnode[which(mapintnode[,1]==gtreel$nodes[n,2]),2]
    }
    if(is.na(gtreel$nodes[gtreel$nodes[n,3],2])){ # internal nodes of a tip in the original tree
      filsd = which(idtips==gtreel$nodes[n,3])
    }else{
      filsd = mapintnode[which(mapintnode[,1]==gtreel$nodes[n,3]),2]
    }
    gtreelO$nodes[idintnodes0,2:3] = c(filsg,filsd)
    gtreelO$nodes[idintnodes0,1] = mapintnode[which(mapintnode[,1]==gtreel$nodes[n,1]),2]
    gtreelO$nodes[filsg,1] = idintnodes0
    gtreelO$nodes[filsd,1] = idintnodes0
    gtreelO$nodes[filsg,4] = gtreel$nodes[gtreel$nodes[n,2],4]
    gtreelO$nodes[filsd,4] = gtreel$nodes[gtreel$nodes[n,3],4]
    idintnodes0 = idintnodes0+1
  }
  # for the root
  filsr = which(gtreelO$nodes[,1]==(ntips+1)) # internal node(s) below root
  if(is.na(gtreel$nodes[gtreel$nodes[rootn,2],2])){ # a tip is directly connected to root in the original tree
    filsg = which(idtips==gtreel$nodes[rootn,2])
    filsd = filsr[1]
    gtreelO$nodes[filsg,1] = ntips+1
    gtreelO$nodes[filsg,4] = gtreel$nodes[gtreel$nodes[rootn,2],4]
    if (length(gtreel$nodes[mapintnode[which(mapintnode[,2]==filsd),1],4])==0){browser()}#ici
    gtreelO$nodes[filsd,4] = gtreel$nodes[mapintnode[which(mapintnode[,2]==filsd),1],4]#ici
  }else if(is.na(gtreel$nodes[gtreel$nodes[rootn,3],2])){
    filsg = filsr[1]
    filsd = which(idtips==gtreel$nodes[rootn,3])
    gtreelO$nodes[filsd,1] = ntips+1
    if (length(gtreel$nodes[mapintnode[which(mapintnode[,2]==filsg),1],4])==0){browser()}#ici
    gtreelO$nodes[filsg,4] = gtreel$nodes[mapintnode[which(mapintnode[,2]==filsg),1],4] ##ici
    gtreelO$nodes[filsd,4] = gtreel$nodes[gtreel$nodes[rootn,3],4]
  }else{
    filsg = filsr[1]
    filsd = filsr[2]
    if(length(which(mapintnode[,2]==filsg))==0 || length(which(mapintnode[,2]==filsd))==0){browser()}
    gtreelO$nodes[filsg,4] = gtreel$nodes[mapintnode[which(mapintnode[,2]==filsg),1],4] #ici
    gtreelO$nodes[filsd,4] = gtreel$nodes[mapintnode[which(mapintnode[,2]==filsd),1],4] #ici
  }
  gtreelO$nodes[ntips+1,2:3] = c(filsg,filsd)
  locationsO = locations[c(idtips,rootn,intnode_ordered),]
  gtreelO$nodes.label = locationsO[,1]
  # reorder so that internal nodes are ordered form oldest to youngest
  # reordered = reorder_treel(gtreelO, locationsO)
  reordered = reorder_treel(gtreelO, locationsO)
  gtreelO = reordered[[1]]
  locationsO = reordered[[2]]

#print(gtreelO)
  check_treel(gtreelO,locationsO)
  return(list(gtreelO,locationsO))
}



# verif function arguments : type, length...
verif_arg <- function(x, type = NA, len = NA, name_function, positive = 0){
	if (sum(is.na(type)) == 0){
		if(!(is.element(class(x),type))){
			stop("in ",name_function ," : invalid 'class' (", class(x), ") for '", deparse(substitute(x)),"'", call. = FALSE)
		} 
	}
	if (sum(is.na(len)) == 0){
		if(!(is.element(length(x),len))){
			stop("in ",name_function ," : parameter '",deparse(substitute(x)),"' is of incorrect size", call. = FALSE)
		}		
	}
	if (positive == 1){
		if(sum(x<0) >0){
			stop("in ",name_function ," : invalid negative value specified for parameter '",deparse(substitute(x)),"'. ","'",deparse(substitute(x)),"' must be positive", call. = FALSE)
		}
	}else if (positive == 2){
		if(sum(x<=0) >0){
			stop("in ",name_function ," : invalid value specified for parameter '",deparse(substitute(x)),"'. ","'",deparse(substitute(x)),"' must be strictly positive", call. = FALSE)
		}	
	}
}



#This function return the minimum value of sigma  wich allow a probability of migration >= at the "proba" parameter
sigma_limit <- function(dmin=20,dmax=500,probalow=.01,probahigh=.99,n_itera=100,sigma_min=0,sigma_max=1000,tries=10000,plotting=0){
  for (i in 1:n_itera){
    testsigma = seq(sigma_min,sigma_max,length.out=tries)
    proba1 = matrix(0,nrow=length(testsigma),ncol=2)
    for(i in 1:length(testsigma)){
      proba1[i,1] = testsigma[i]
      proba1[i,2]=dnorm(dmin,0,sqrt(testsigma[i]))/(dnorm(0,0,sqrt(testsigma[i]))) 
    }
  }
  for (i in 2:tries) {
    sig_lower_limit=proba1[i,1]
    if(proba1[i,2]>=probalow){
      break
    }
  }
  for (i in 1:n_itera){
    proba2 = matrix(0,nrow=length(testsigma),ncol=2)
    testsigma2=seq(sig_lower_limit,1000000,length.out=tries)
    for(i in 1:length(testsigma2)){
      proba2[i,1] = testsigma2[i]
      proba2[i,2]=dnorm(dmax,0,sqrt(testsigma2[i]))/(dnorm(0,0,sqrt(testsigma2[i])))
    }
  }
  for (i in 2:tries) {
    sig_upper_limit=proba2[i,1]
    if(proba2[i,2]>0.9){
      break
    }
  }
  return(c(sig_lower_limit,sig_upper_limit))
}

sigma_limit_V2 <- function(dist,proba_lim = 0.1){
  start=0
  Borne=1e25
  testsigma2=seq(start,Borne,length.out = 100)
  j=1
  for (x in 1:100) {
    while(j<=100)
    {
      proba2=dnorm(dist,0,sqrt(testsigma2[j]))/(dnorm(0,0,sqrt(testsigma2[j])))
      if(proba2>=proba_lim){
        sig_limit=testsigma2[j]
        testsigma2=seq(testsigma2[j-1],testsigma2[j],length.out = 100)
        j=1
        break
      }
      else{
        j=j+1
      }
    }
  }
  return(sig_limit)
}

