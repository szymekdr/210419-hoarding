# loading required packages ----
library(pacman)
p_load(brms, phytools, geiger, tidyverse, ape)

library(here)

options(mc.cores = parallel::detectCores())

# load data and tree ----
data <- read_csv(here('data', 'hoarding_data_27Jul2020.csv'))
tree <- read.tree(here('data', 'Stage2_Hackett_MCC_no_neg.tre'))
tree_final <- drop.tip(tree, setdiff(tree$tip.label, data$tip_label))

data_analysis <- filter(data, hoarding_status != 'unknown_hoarding_status')
data_analysis <- as.data.frame(data_analysis)

names(data_analysis)[which(names(data_analysis) == 'family')] <- 'tax_family'


# define brms model ----
tree_final$edge.length[tree_final$edge.length == 0] <- 1e-10
A <- vcv.phylo(tree_final)


model <- brmsformula(
  hoarding_status ~ 1 + (1|gr(tip_label, cov = A))
)

# get_prior(model, data = data_analysis, family = categorical, data2 = list(A = A))

modelrun <- brm(model, chains = 2, iter = 10000, warmup = 3000, thin = 5,
                data = data_analysis,
                data2 = list(A = A),
                family = categorical,
                control = list(adapt_delta = 0.9, max_treedepth = 15)
                )

