# loading required packages ----
library(pacman)
p_load(MCMCglmm, phytools, geiger, tidyverse, ape)


# load data and tree ----
data <- read_csv(here('data', 'hoarding_data_27Jul2020.csv'))
tree <- read.tree(here('data', 'Stage2_Hackett_MCC_no_neg.tre'))

data_analysis <- filter(data, hoarding_status != 'unknown_hoarding_status')
data_analysis <- as.data.frame(data_analysis)

names(data_analysis)[which(names(data_analysis) == 'family')] <- 'tax_family'


# define and run the mcmcglmm model ----
j <- length(unique(data_analysis$hoarding_status))
IJ <- (1/j)*(diag(j-1) + matrix(1, j-1, j-1))

# making sure all edges are non zero but tree remians ultrameric
tree_final$edge.length[tree_final$edge.length == 0] <- 1e-10

# form phylogenetic VCV scaled to unity r and produce (sparse) inverse of A
phylo_vcv <- inverseA(tree_final, nodes = 'ALL')

prior2 <- list(B = list(mu = rep(0, j-1), V = kronecker(IJ, 1) * (1.7 +pi^2/3)),
               R = list(V = IJ, fix = 1),
               G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4)))
model2 <- MCMCglmm(hoarding_status ~ trait - 1, 
                   random = ~tip_label,
                   rcov = ~us(trait):units,
                   prior = prior2,
                   ginverse = list(tip_label = phylo_vcv$Ainv),
                   data = data_analysis, family = 'categorical',
                   verbose = F,
                   nitt = 1e7, burnin = 1e6, thin = 5000)


# collect and save relevant objects ----
save(model2, file = 'out_model2_mcmcglmm.Rdata')