#### Correct Model

library(lavaan)  
library(Amelia)
library(mice)
library(semTools)
library(methods) 
library(gtools)  

################## Conditions #####################
Meth <- "CorrectModel" 

NonNorm <- c(1, 2, 3)   ### 1 = (S = 1.5, K = 3); 2 = (S = 2, K = 7); 3 = (S = 3, K = 21)
N <- c(150, 300, 600)   ### c(150, 300, 600)  
MissMec <- c(1, 2, 3)   ### 1 = MCAR; 2 = MAR-Head; 3 = MAR-Tail
MissPro <- c(0.15, 0.3) ### c(0.15, 0.3)

ImpMeth <- c(1, 2) ### 1 = Amelia, 2 = mice::pmm
EstMeth <- c(1, 2, 3) ### 1 = MLR, 2 = MLM, 3 = MLMV

cond.list <- expand.grid(Meth = Meth, NonNormn = NonNorm, N = N, MissMec = MissMec, MissPro = MissPro, 
						ImpMeth = ImpMeth, EstMeth = EstMeth)

#SK <- list(c(1.5, 3), c(2, 7), c(3, 21))
###################################################
set.seed(980210)
seedList <- sample(1:999999, 1000)

#####~~~~~~~~~~~test~~~~~~~~~~~~~~~#####
# cond <- cond.list[1,]
# SUBNUM <- 4

# cond.list <- cond.list[1:2,]
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

nRep <- 1000

for (SUBNUM in 1:nRep) {
RMIFUN <- function (cond){ 
	Meth <- as.integer(cond[1])
	NonNorm <- as.integer(cond[2])
	N <- as.integer(cond[3]) 
	MissMec <- as.integer(cond[4])
	MissPro <- as.numeric(cond[5]) 
	ImpMeth <- as.numeric(cond[6]) 
	EstMeth <- as.numeric(cond[7]) 
	
	RESNAME_POP <- paste0("nonnorm", NonNorm) 
	
	pop.nonnorm <- read.csv(paste0(getwd(),"/Population/", RESNAME_POP, ".csv"), header = TRUE) 
	#head(pop.nonnorm)
	
	model <- "
	  #### measurement model
		f1 =~ x1 + x2 + x3
		f2 =~ x4 + x5 + x6
		f3 =~ x7 + x8 + x9
	  #### regressions
		f2 ~ f1
		f3 ~ f1 + f2
		
		f1 ~~ f1
		f2 ~~ f2
		f3 ~~ f3
	"
	#model.pop <- sem(model, data = pop.nonnorm, estimator = "MLMV")
	#summary(model.pop)
	
  #############################################################################################
  #### Missing data
	
	######## patterns 
	pattern1 <- combinations(6, 1, v = c(1,2, 4, 5, 7, 8))
	pattern2 <- combinations(6, 2, v = c(1,2, 4, 5, 7, 8))
	pattern3 <- combinations(6, 3, v = c(1,2, 4, 5, 7, 8))
	pattern4 <- combinations(6, 4, v = c(1,2, 4, 5, 7, 8))
	pattern5 <- combinations(6, 5, v = c(1,2, 4, 5, 7, 8))
	pattern6 <- combinations(6, 6, v = c(1,2, 4, 5, 7, 8))
	
	makepatternfun <- function(pat, ncolumn){
		patMatrix <- matrix(1, nrow = nrow(pat), ncol=ncolumn)

		for (i in 1:nrow(pat)){
			patMatrix[i, pat[i,]] <- 0 
		}
		return(patMatrix)
	}
	patMatrix1 <- makepatternfun(pattern1, 9)
	patMatrix2 <- makepatternfun(pattern2, 9)
	patMatrix3 <- makepatternfun(pattern3, 9)
	patMatrix4 <- makepatternfun(pattern4, 9)
	patMatrix5 <- makepatternfun(pattern5, 9)
	patMatrix6 <- makepatternfun(pattern6, 9)
	
	patMatrix <- rbind(patMatrix1, patMatrix2, patMatrix3,
					 patMatrix4, patMatrix5, patMatrix6)
	
	######## weights
	myWeights <- patMatrix
	myWeights[,] <- 0

	myWeights[patMatrix[,1]==0 | patMatrix[,2]==0, 3] <- 1
	myWeights[patMatrix[,4]==0 | patMatrix[,5]==0, 6] <- 1
	myWeights[patMatrix[,7]==0 | patMatrix[,8]==0, 9] <- 1
	
	######## type
	myType2 <- rep("RIGHT", nrow(patMatrix))
	myType3 <- rep("LEFT", nrow(patMatrix))

	######## amputation
	if (MissMec==1) {
	  pop.nonnorm.miss <- ampute(data = pop.nonnorm, prop = MissPro*2, patterns = patMatrix, mech = "MCAR")$amp
	} else if (MissMec==2) {
	  pop.nonnorm.miss <- ampute(data = pop.nonnorm, prop = MissPro*2, patterns = patMatrix, mech = "MAR", weights = myWeights, type = myType2)$amp
	} else if (MissMec==3) {
	  pop.nonnorm.miss <- ampute(data = pop.nonnorm, prop = MissPro*2, patterns = patMatrix, mech = "MAR", weights = myWeights, type = myType3)$amp
	}
	####### MissPro*2, because here prop = the proportion of incomplete rows, not the same as the proportion of missing cells
	
	#mice::md.pattern(pop.nonnorm.miss, plot = FALSE)

	#############################################################################################
	#### Sampling
	set.seed(seedList[SUBNUM]) 
	cont.miss <- pop.nonnorm.miss[sample(1:nrow(pop.nonnorm.miss), N),]
	
	#mice::md.pattern(cont.miss, plot = FALSE)

	#############################################################################################
	#### Impute and Combine
  source("2.RobustFitMI_functions.R")

	NImp <- 20

	out <- MethAll(dat = cont.miss, meth = Meth, impseed = seedList[SUBNUM], 
	                  nimp = NImp, impmeth = ImpMeth, estmeth = EstMeth)
	
	return(out)
}

result <- apply(cond.list, 1, RMIFUN)
dput(result, paste0("RobustFItMI_meth", Meth, "_rep.", SUBNUM, ".Rdata"))
}
