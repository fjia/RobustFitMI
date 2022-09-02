#### Foldnes and Olsson (2016) mehtod
	
######### marker variable
library(lavaan)  

latent.b1 <- matrix(c(0.4, 0.286), 2, byrow = TRUE)  ##x

latent.b23 <- matrix(c(0, 0, 0.286, 0), 2, byrow = TRUE)  ##y
							
latent.rsv <- matrix(c(0.412, 0, 0, 0.378), 2)


latent.sigma.yy <- solve(diag(2) - latent.b23) %*% (0.49 * latent.b1 %*% t(latent.b1) + latent.rsv) %*% t(solve(diag(2) - latent.b23))
latent.sigma.xx <- 0.49
latent.sigma.yx <- solve(diag(2) - latent.b23) %*% (0.49 * latent.b1)
latent.sigma.xy <- 0.49 * t(latent.b1) %*% t(solve(diag(2) - latent.b23))

latent.sigma.left <- rbind(latent.sigma.yy, latent.sigma.xy)
latent.sigma.right <- rbind(latent.sigma.yx, latent.sigma.xx)

latent.sigma <- cbind(latent.sigma.left, latent.sigma.right)

colnames(latent.sigma) <- c("f2", "f3", "f1")
rownames(latent.sigma) <- c("f2", "f3", "f1")
lambda.mx <- matrix(0, 9, 3) 
lambda.mx[1:3, 1] <- 1  ##0.7
lambda.mx[4:6, 2] <- 1  ##0.7
lambda.mx[7:9, 3] <- 1  ##0.7

res.var <- rep(0.51, 9)

obs.sigma <- lambda.mx %*% latent.sigma %*% t(lambda.mx)  + diag(res.var)  ### pop cov
colnames(obs.sigma) <- c("x4", "x5", "x6","x7", "x8", "x9","x1", "x2", "x3")
obs.rho <- cov2cor(obs.sigma)  ### pop cor

npop <- 500000
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
# fit <- sem(model, sample.cov = obs.sigma, sample.nobs = npop)
# summary(fit, standardized = TRUE)

source("gen.nonnormal.fo2016.r") #### Foldnes and Olsson (2016)
SK <- list(c(1.5, 3), c(2, 7), c(3, 21))
NonNorm <- c(1, 2, 3)  ### 1 = (S = 1.5, K = 3), 2 = (S = 2, K = 7), ### 3 = (S = 3, K = 21).
cond.list.POP <- expand.grid(NonNorm = NonNorm)

dir.create(paste0(getwd(),"/Population/"))

#####~~~~~~~~~~~test~~~~~~~~~~~~~~~#####
# cond <- cond.list.POP[1,]
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

nonnormPopFUN <- function (cond){ 
  NonNorm <- as.integer(cond[1])
	RESNAME_POP <- paste0("nonnorm", NonNorm) 
	
	sk <- rep(SK[[NonNorm]][1], 9) #desired level of skewness for each variable
	ku <- rep(SK[[NonNorm]][2], 9) # desired level of kurtosis for each variables
	set.seed(980210)
	pop.nonnorm <- gen.nonnormal.f02016(npop, obs.sigma, sk, ku, 9) 
	colnames(pop.nonnorm) <- colnames(obs.sigma)
	write.csv(pop.nonnorm, file=paste0(getwd(), "/Population/", RESNAME_POP, ".csv"), row.names=FALSE, quote = FALSE)
}

apply(cond.list.POP, 1, FUN =nonnormPopFUN)


