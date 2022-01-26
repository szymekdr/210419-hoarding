library(phylolm)
#library(phytools)
#library(mvMORPH)
library(shiftPlot)
library(MuMIn)
library(brms)

setwd("~/Documents/Work/Research/HoardingBirds")

#there have been new modifications to the file I saved out in previous scripts. load
#it in here
dat <- read.csv("data/hoarding_data_27Jul2020.csv", stringsAsFactors=FALSE)

#read in the big bird tree
jetz <- read.tree("data/Stage2_Hackett_MCC_no_neg.tre")

#create a coloniality column
dat$colonial <- 1
dat$colonial[dat$INTRA_COL=="acolonial"] <- 0

#create a family/cooperative or no family breeding system column
dat$family.breeding <- 1
dat$family.breeding[dat$fam_sys_inferred_nk50=="no_fam"] <- 0

#come up with a carnivore column that is carnivore+vertivore
dat$our.carnivore <- dat$carnivore+dat$vertivore

#set row names for phylolm model
row.names(dat) <- dat$tip_label

#cut the tree down to spp in dat
pruned <- ladderize(drop.tip(jetz, tip=setdiff(jetz$tip.label, dat$tip_label)))

#let's try a model. here's the original, first big result we shared. 
OUmodel <- phylolm(hoard~colonial+invertivore+aquatic+our.carnivore+nectarivore+frugivore+granivore+Bio1_mean+Bio12_mean,
	data=dat, phy=pruned, model="OUrandomRoot")

#here's one that mario michael suggested more recently. prune down
datCut <- dat[,c("hoard","colonial","family.breeding","median_clutch_size","mov_max","body_weight",
	"food_energy_mean","foodPCA1","foodPCA2","foodPCA3","foodPCA4","Bio1_mean","Bio12_mean","strat_final")]
datCut <- datCut[complete.cases(datCut),]

#log two of the variables
datCut$log.mass <- log(datCut$body_weight)
datCut$log.clutch <- log(datCut$median_clutch_size)

#scale and center the variables
setAside <- scale(datCut[,c("log.mass","log.clutch","mov_max",
	"food_energy_mean", "foodPCA1","foodPCA2","foodPCA3", "foodPCA4","Bio1_mean","Bio12_mean")],
	scale=TRUE, center=TRUE)
finalDat <- data.frame(setAside, hoard=datCut$hoard, hoarding.strategy=datCut$strat_final,
	family.breeding=datCut$family.breeding, colonial=datCut$colonial)

#prune tree
pruned2 <- drop.tip(pruned, tip=setdiff(pruned$tip.label, row.names(finalDat)))

#run some models
OUmodel2 <- phylolm(hoard~colonial+log.clutch+mov_max+log.mass+
	food_energy_mean+foodPCA1+foodPCA2+foodPCA3+foodPCA4+Bio1_mean+Bio12_mean,
	data=finalDat, phy=pruned2, model="OUrandomRoot")


######SKIP THIS


################################################################################
########################### TRY SOME BRMS MODELS ###############################
################################################################################

#create a covariance matrix
A <- vcv(pruned2)

#create ordered factors for the model
forBRMS <- finalDat
forBRMS$species <- row.names(forBRMS)
<- factor(forBRMS$hoarding.strategy, ordered=TRUE)

model_simple <- brm(the.factor ~ Bio1_mean + (1|gr(species, cov = A)),
	data = forBRMS, 
	family = cumulative("logit"),
	data2 = list(A = A),
	prior = c(
	set_prior("normal(0,3)", class = "b"),
	set_prior("normal(0,3)", class = "Intercept"),
	set_prior("normal(0,4)", class = "sd", coef = "Intercept", group = "species"))
)






#split on hoarding type
prov <- finalDat[finalDat$hoarding.strategy=="provision" | finalDat$hoarding.strategy=="no",]
prunedProv <- drop.tip(pruned2, tip=setdiff(pruned2$tip.label, row.names(prov)))

sust <- finalDat[finalDat$hoarding.strategy=="sustenance" | finalDat$hoarding.strategy=="no",]
prunedSust <- drop.tip(pruned2, setdiff(pruned2$tip.label, row.names(sust)))

provModel <- phylolm(hoard~family.breeding+colonial+log.clutch+mov_max+log.mass+
	food_energy_mean+foodPCA1+foodPCA2+foodPCA3+foodPCA4+Bio1_mean+Bio12_mean,
	data=prov, phy=prunedProv, model="lambda")
sustModel <- phylolm(hoard~colonial+log.clutch+mov_max+log.mass+
	food_energy_mean+foodPCA1+foodPCA2+foodPCA3+foodPCA4+Bio1_mean+Bio12_mean,
	data=sust, phy=prunedSust, model="lambda")

#provDredge <- dredge(phylolm(hoard~colonial+family.breeding+log.clutch+mov_max+log.mass+
#	food_energy_mean+foodPCA1+foodPCA2+foodPCA3+foodPCA4+Bio1_mean+Bio12_mean,
#	data=prov, phy=prunedProv, model="lambda"), trace=2)
#saveRDS(provDredge, "data/provDredge.RDS")
provDredge <- readRDS("data/provDredge.RDS")

#sustDredge <- dredge(phylolm(hoard~colonial+family.breeding+log.clutch+mov_max+log.mass+
#	food_energy_mean+foodPCA1+foodPCA2+foodPCA3+foodPCA4+Bio1_mean+Bio12_mean,
#	data=sust, phy=prunedSust, model="lambda"), trace=2)
#saveRDS(sustDredge, "data/sustDredge")
sustDredge <- readRDS("data/sustDredge.RDS")

provModav <- model.avg(provDredge, subset = delta < 2)
sustModav <- model.avg(sustDredge, subset = delta < 2)





#####START AGAIN HERE





#corHMM tests. sort data
prepped <- dat[pruned$tip.label,c("tip_label","hoard")]

#now just dive right into the vignette code
#(https://rdrr.io/github/thej022214/corHMM/f/vignettes/corHMMv2.0-vignette.Rmd)
Precur_LegendAndMat <- getStateMat4Dat(prepped)
Precur_LegendAndMat
Precur_R1 <- Precur_LegendAndMat$rate.mat
Precur_R1 <- dropStateMatPars(Precur_R1, c(1,2))
Precur_R1
Precur_R2 <- Precur_LegendAndMat$rate.mat
#I think if you hash out the line below it'll run an all rates different model
Precur_R2 <- equateStateMatPars(Precur_R2, c(1,2))
Precur_R2
RateClassMat <- getRateCatMat(2) #
RateClassMat <- equateStateMatPars(RateClassMat, c(1,2))
RateClassMat
Precur_FullMat <- getFullMat(list(Precur_R1, Precur_R2), RateClassMat)
Precur_FullMat[c(4,2), c(2,4)] <- 0
Precur_FullMat
plotMKmodel(Precur_FullMat, 2, display = "row", text.scale = 0.7)
Precur_res.corHMM <- corHMM(phy = pruned, data = prepped, rate.cat = 2, rate.mat = Precur_FullMat, nstarts=3)
Precur_res.corHMM
plotMKmodel(Precur_res.corHMM, display = "row", text.scale = 0.7)

#now the more complicated task of plotting the pie charts. make a simple function
nodePlotter <- function(tree, nodeVals, interest.col, threshold, small.size, large.size, cols)
{
	#plotrix has this fairly insane behavior where if one of the numbers in the vector is 0,
	#it will simply skip it, which means circles of 100% probability get all fucked up. add a
	#teeny fudge factor so this doesn't happen
	nodeVals[,1][nodeVals[,1]==0] <- 0.001
	nodeVals[,2][nodeVals[,2]==0] <- 0.001

	#for(i in 1:1)
	for(i in 1:dim(nodeVals)[1])
	{
		#have a running node counter so i varies with the actual node numbers
		nodeNumbers <- row.names(nodeVals)
		node <- as.numeric(nodeNumbers[i])

		#check out the node. if the column of interest is > than the threshold,
		#note it here. this is counter intuitive, but this pulls a single row which automatically
		#gets converted to a vector, so there is no comma below before interest.col
		temp <- nodeVals[i,]

		if(temp[interest.col] > threshold)
		{
			pieSize <- large.size
		}
		else
		{
			pieSize <- small.size
		}

		#lay the pie chart down
		floating.pie(xpos=node.depth.edgelength(tree)[node], ypos=node.height(tree)[node],
			x=temp, radius=pieSize, col=cols, border=NA)
	}
}

#set the nodes up to look right for this function. the sum of the probability that the node
#is in either state of nectarivory, but likely in rate 2, seems like what you want

###### ONE OF THE FOLLOWING OBJECTS IS PROBABLY WHAT YOU WANT FOR ASRs. 
#maybe this next obj right there. "ones" refers to presence of precursor



pulledNodes <- data.frame(zero=1, one=Precur_res.corHMM$states[,3]+Precur_res.corHMM$states[,4])
pulledNodes$zero <- pulledNodes$zero-pulledNodes$one
pulledNodes <- as.matrix(pulledNodes)
colnames(pulledNodes) <- c(0,1)
rownames(pulledNodes) <- (length(pruned$tip.label)+1):(length(pruned$tip.label)*2-1)

tipColors <- rep("black", length(pruned$tip.label))
tipColors[prepped$hoard=="1"] <- "red"
plot(pruned, cex=0.1, tip.col=tipColors)
nodePlotter(tree=pruned, nodeVals=pulledNodes, interest.col=2, threshold=0.05, small.size=0.1,
	large.size=0.9, cols=c("black","red"))

cexSize <- rep(0.01,length(tipColors))
cexSize[tipColors=="red"] <- 0.05

pdf(file="outputs/ASR.pdf", width=10, height=10)
plot(pruned, cex=cexSize, tip.col=tipColors, edge.width=0.05)
nodePlotter(tree=pruned, nodeVals=pulledNodes, interest.col=2, threshold=0.05, small.size=0.05,
	large.size=0.2, cols=c("black","red"))
dev.off()




