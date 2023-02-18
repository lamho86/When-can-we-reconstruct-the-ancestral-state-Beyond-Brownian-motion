# Check for required packages
for(pack in c("ape", "geiger", "phytools")){
    if(pack %in% rownames(installed.packages()) == F){
        install.packages(pack)
    }
}

library(ape)
library(geiger)
library(phytools)
# Source code of Sim_BBM.R & BBM_functions_bounds_estimated_or_not_plus_uncertainty_with_CIs.R
# ... can be clone from this repo: http://github.com/fcboucher/BBM
source("Sim_BBM.R")
source("BBM_functions_bounds_estimated_or_not_plus_uncertainty_with_CIs.R")

generate_forest <- function(N, seq_gen, rseed = 1){
    forest <- list()
    set.seed(rseed)
    tree <- rtree(2, rooted = TRUE)
    tree$edge.length <- seq_gen(1) * rep(1, length(tree$edge.length))
    if(2 %in% N){
        forest[["2"]] <- tree
    }
    
    i <- 2
    while(i <= max(N)){
        label = paste("t", i + 1, sep = "")
        tree <- bind.tip(
            tree, tip.label = label,
            where = i, position = seq_gen(i)
        )
        
        if((i + 1) %in% N){
            forest[[paste(i + 1)]] <- tree
        }
        
        i <- i + 1
    }
    return(forest)
}

estimation <- data.frame() # Empty data frame

leaves <- c(10, 50, 100, 500, 1000, 5000)
forest <- generate_forest(leaves, function(n){return(1/n)})

for(i in 1:1000){
    for(j in leaves){
        tree <- forest[[paste(j)]]
        trait = Sim_BBM(
            tree, x0 = 0, Npts = 50, 
            sigma = sqrt(1/2), bounds = c(-1, 1)
        )
        BBM = fit_BBM_model_uncertainty(
            tree,
            trait = trait, Npts = 50, 
            bounds = c(-1, 1), uncertainty = F
        )
        estimation <- rbind(estimation, c(j, BBM$par$root_value))
    }
}
names(estimation) <- c("num_leaves", "estimates")
write.csv(estimation, "BBM_notBigBang.csv")
