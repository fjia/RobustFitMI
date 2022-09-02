

MethAll <- function(dat, meth, impseed, nimp, impmeth, estmeth){
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
	if (impmeth==1){
		#### Amelia imputation
		imp <- Amelia::amelia(dat, m=nimp, p2s=0, empri =.01*nrow(dat), seed = impseed)
		dat.imp <- list()
		for (i in 1:nimp){
			dat.imp[[i]] <- imp$imputations[[i]]
		}
	} else {			
		#### mice-PMM imputation
		imp <-  mice::mice(dat, meth=c(rep("pmm",2),"",rep("pmm",2),"",rep("pmm",2),""), m=nimp ,print=FALSE, maxit=20, seed = impseed)
		imputations <- mice::complete(imp, "long")
		dat.imp <- list()
		for (i in 1:nimp){
			dat.imp[[i]] <- imputations[imputations$.imp==i,-1:-2]
			dat.imp[[i]] <- data.frame(do.call("cbind", dat.imp[[i]]))
		}
	}
	
	if (estmeth==1){
		fit.list <- semTools::sem.mi(model, data = dat.imp, estimator = "MLR")
	} else if (estmeth==2){
		fit.list <- semTools::sem.mi(model, data = dat.imp, estimator = "MLM")
	} else {
		fit.list <- semTools::sem.mi(model, data = dat.imp, estimator = "MLMV")	
	}
	
	#### Pool Test Stat and fit
	test.D3 <- try(semTools::lavTestLRT.mi(fit.list, test = "D3", pool.robust = FALSE, asymptotic = TRUE)) ##D3SN	
	test.D2 <- try(semTools::lavTestLRT.mi(fit.list, test = "D2", pool.robust = FALSE, asymptotic = TRUE)) ##D2ASN
	test.D2R <- try(semTools::lavTestLRT.mi(fit.list, test = "D2", pool.robust = TRUE)) ##D2
	test.D2RA <-try(semTools::lavTestLRT.mi(fit.list, test = "D2", pool.robust = TRUE, asymptotic = TRUE)) ##D2A
	
	#####################################
	### 1)
	test.d3.chi <- test.D3["chisq.scaled"]   
	test.d3.df <- test.D3["df.scaled"] 
	test.d3.p <- test.D3["pvalue.scaled"]
	
	### 2)
	test.d2.chi <- test.D2["chisq.scaled"]   
	test.d2.df <- test.D2["df.scaled"] 
	test.d2.p <- test.D2["pvalue.scaled" ] 
	
	### 3)
	test.d2r.f <- test.D2R["F.scaled"]   
	test.d2r.df1 <- test.D2R["df1.scaled"]
	test.d2r.df2 <- test.D2R["df2.scaled"]  
	test.d2r.p <- test.D2R["pvalue.scaled"] 
	
	### 4)
	test.d2ra.chi <- test.D2RA["chisq.scaled"]   
	test.d2ra.df <- test.D2RA["df.scaled"] 
	test.d2ra.p <- test.D2RA["pvalue.scaled"] 
	

	final <- c(test.d3.chi=test.d3.chi, test.d3.df=test.d3.df, test.d3.p=test.d3.p,
	           test.d2.chi=test.d2.chi, test.d2.df=test.d2.df, test.d2.p=test.d2.p,
	           test.d2r.f=test.d2r.f, test.d2r.df1=test.d2r.df1, test.d2r.df2=test.d2r.df2, test.d2r.p=test.d2r.p,
	           test.d2ra.chi=test.d2ra.chi, test.d2ra.df=test.d2ra.df, test.d2ra.p=test.d2ra.p)
	
	return(final) 
}