#
#
# Bayesian zero-truncated negative binomial generalized linear mixed model for count data
#
# Function name: constraints.mcmc
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 22 AUG 2016
#
# Model statement:
#	z_ij|z_ij>0 ~ ZTNB(lambda_ij,alpha)
#	log(lambda_ij) = x_ij%*%beta_j
# 	beta_j ~ N(mu_beta,Sigma)
#	mu_beta ~ N(0,sigma.beta^2*I)
#	Sigma ~ Wish(S_0,nu)
#	alpha ~ Gamma(a,b)
#	note: E[z_ij]=lambda_ij and Var[z_ij]=lambda_ij+lambda_ij^2/alpha
#	note: E[alpha]=a/b and Var[alpha]=a/(b^2)
#
# Reference:
#
# Required R packages: mvtnorm (if using multivariate normal prior on beta)
#
# Inputs:
#
# z - Vector of length n containing the count of observations corresponding to 
# 	each row in the design matrix X. Note that the value of z[1] corresponds to
#	X[1,], z[2] corresponds to X[2,], etc. Also note that z>0.
# X - Design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# g - variable that defines groups of observations in z
# priors - list of priors containing the following elements:
#	1. sigma.beta - Standard deviation of normal prior on mu.beta
#	2. S0 - Scale matrix for the inverse-Wishart prior on Sigma
#	3. nu - Degrees of freedom for the IW prior on Sigma
#	4. a - Shape parameter of gamma prior for alpha 
#	5. b - Rate parameter of gamma prior for alpha
# start - list of starting values containing the following elements:
#	1. beta - Vector of starting values for coefficients
#	2. mu.beta - Vector of starting values for mean of betas
#	3. Sigma - Variance-covariance matrix for betas
#	4. alpha - Over-dispersion parameter for negative binomial distribution
# tune - List of tuning parameters containing the following elements:
#	1. beta - Tuning parameter for Metropolis-Hastings update on beta
#	2. alpha - Tuning parameter for Metropolis-Hastings update on alpha
# adapt - Switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - Number of desired MCMC iterations
#
#

constraints.mcmc <- function(s,S,priors,start,tune,n.mcmc,adapt=TRUE){
		
	t.start <- Sys.time()
	cat(paste("Start time:",t.start,"\n"))
		
	
	###
	### Libraries and Subroutines
	###
	
	library(mvtnorm) #For beta propsals
	library(raster) #For spatial operations (e.g., using S as a mask for mu proposals)
	
	get.dmvt2 <- function(x,y,sigma,a,r,nu){  # calculate log density of mixture t-distribution
		d <- 2
		E.1 <- sigma^2*matrix(c(1,sqrt(a)*r,sqrt(a)*r,a),2,2,byrow=TRUE)
		E.2 <- sigma^2*matrix(c(1,sqrt(a)*-r,sqrt(a)*-r,a),2,2,byrow=TRUE)	
		c <- lgamma((nu+d)/2)-(lgamma(nu/2)+log(nu)+log(pi))+log(det(E.1)^-0.5)
		dif <- x-y
		# try using mahalanobis() instead?
		kernel.1 <- (1+1/nu*(rowSums((dif%*%solve(E.1))*dif)))^-((nu+d)/2)
		kernel.2 <- (1+1/nu*(rowSums((dif%*%solve(E.2))*dif)))^-((nu+d)/2)
		log(0.5)+c+log((kernel.1+kernel.2))
	}

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}


	###
	###  Setup Variables 
	###
	
	# browser()
	T <- nrow(s)  # number of telemetry locations	
	
	# Index of locations in each error class
	
	lc <- as.numeric(s$lc)  # error classes

	lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)
	s.list <- lapply(lc.list,function(x) s@coords[x,])
	T.list <- lapply(lc.list,length)
	m <- length(lc.list)  # number of error classes

	mu.star <- matrix(0,T,2)
		
	s <- s@coords  # coordinates of telemetry locations
				
	###
	### Starting values
	###
	
	mu <- start$mu  # true animal locations
	# cell.idx <- extract(S,mu,cellnumbers=TRUE)[,1]  # cells in S corresponding to mu

	# observation model parameters
	sigma <- start$sigma
	a <- start$a
	rho <- start$rho
	nu <- start$nu


	###
	### Standardize parameters
	###

	# cat("\nStandardizing variables....")
	# s.scale <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# s.center <- apply(s,2,function(x) max(x)+min(x))/2

	# s <- (s-matrix(s.center,nrow=N,ncol=2,byrow=TRUE))/s.scale
	# mu <- (mu-matrix(s.center,nrow=N,ncol=2,byrow=TRUE))/s.scale
	# sigma <- sigma/s.scale
	
		
	###
	### Priors and tuning parameters
	###

	# priors$sigma <- priors$sigma/s.scale
	# tune$mu <- tune$mu/s.scale
	# tune$sigma <- tune$sigma/s.scale


	###
	### Receptacles for output
	###

	cat("\nCreating receptacles for output....")
	mu.save <- array(0,dim=c(T,2,n.mcmc))
	sigma.save <- matrix(0,n.mcmc,m)
	a.save <- matrix(0,n.mcmc,m)
	rho.save <- matrix(0,n.mcmc,m)
	nu.save <- matrix(0,n.mcmc,m)
	
	keep <- list(sigma=rep(0,m),a=rep(0,m),rho=rep(0,m),nu=rep(0,m),mu=rep(0,3))
	keep.tmp <- keep  # for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning
	
	###
	###  Begin MCMC Loop 
	###

	t.mu.update <- 0
	t.mcmc.start <- Sys.time()
	cat("\nEntering MCMC Loop....\n")

	for(k in 1:n.mcmc){ #Loop to iterate MCMC algorithm: Appendix B, step 7
		if(k%%100==0) cat(k,"");flush.console()

		if(adapt==TRUE & k%%Tb==0) {  # adaptive tuning
# browser()
			keep.tmp$mu <- keep.tmp$mu/(unlist(T.list)*Tb)
			keep.tmp$sigma <- keep.tmp$sigma/Tb
			keep.tmp$a <- keep.tmp$a/Tb
			keep.tmp$rho <- keep.tmp$rho/Tb
			keep.tmp$nu <- keep.tmp$nu/Tb
			
			tune$sigma <- sapply(1:m,function(x) get.tune(tune$sigma[x],keep.tmp$sigma[x],k))
			tune$a <- sapply(1:m,function(x) get.tune(tune$a[x],keep.tmp$a[x],k))
			tune$rho <- sapply(1:m,function(x) get.tune(tune$rho[x],keep.tmp$rho[x],k))
			tune$nu <- sapply(1:m,function(x) get.tune(tune$nu[x],keep.tmp$nu[x],k))
			tune$mu <- sapply(1:m,function(x) get.tune(tune$mu[x],keep.tmp$mu[x],k))
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

			
		###
		###  Sample components of Sigma
		###

		sigma.star <- rnorm(m,sigma,tune$sigma)  # proposals for sigma
		a.star <- rnorm(m,a,tune$a)  # proposals for a
		rho.star <- rnorm(m,rho,tune$rho)  # proposals for rho
		nu.star <- rnorm(m,nu,tune$nu)  # proposals for nu
		
# browser()				
		for(i in 1:m){ #Loop to iterate over error classes: Appendix B, step 2(f)
# i <- 1
			lc.idx <- lc.list[[i]] #Index of locations in error class i
			s.tmp <- s[lc.idx,]
			mu.tmp <- mu[lc.idx,]
			T.tmp <- T.list[[i]]
				
			### Sample sigma: Appendix B, step 2(b)

			if(sigma.star[i]>priors$sigma[1] & sigma.star[i]<priors$sigma[2]){
				mh.star.sigma <- sum(get.dmvt2(s.tmp,mu.tmp,sigma.star[i],a[i],rho[i],nu[i]))
				mh.0.sigma <- sum(get.dmvt2(s.tmp,mu.tmp,sigma[i],a[i],rho[i],nu[i]))
				if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					sigma[i] <- sigma.star[i]
					keep$sigma[i] <- keep$sigma[i]+1
					keep.tmp$sigma[i] <- keep.tmp$sigma[i]+1
				}
			}

			### Sample a

			if(a.star[i]>priors$a[1] & a.star[i]<priors$a[2]){
				mh.star.a <- sum(get.dmvt2(s.tmp,mu.tmp,sigma[i],a.star[i],rho[i],nu[i]))
			  	mh.0.a <- sum(get.dmvt2(s.tmp,mu.tmp,sigma[i],a[i],rho[i],nu[i]))
				if(exp(mh.star.a-mh.0.a)>runif(1)){
					a[i] <- a.star[i]
					keep$a[i] <- keep$a[i]+1
					keep.tmp$a[i] <- keep.tmp$a[i]+1
				}
			}

			### Sample rho

			if(rho.star[i]>priors$rho[1] & rho.star[i]<priors$rho[2]){
				mh.star.rho <- sum(get.dmvt2(s.tmp,mu.tmp,sigma[i],a[i],rho.star[i],nu[i]))
			  	mh.0.rho <- sum(get.dmvt2(s.tmp,mu.tmp,sigma[i],a[i],rho[i],nu[i]))
				if(exp(mh.star.rho-mh.0.rho)>runif(1)){
					rho[i] <- rho.star[i]
					keep$rho[i] <- keep$rho[i]+1
					keep.tmp$rho[i] <- keep.tmp$rho[i]+1
				}
			}
			
			### Sample nu

			if(nu.star[i]>priors$nu[1] & nu.star[i]<priors$nu[2]){
				mh.star.nu <- sum(get.dmvt2(s.tmp,mu.tmp,sigma[i],a[i],rho[i],nu.star[i]))
		  		mh.0.nu <- sum(get.dmvt2(s.tmp,mu.tmp,sigma[i],a[i],rho[i],nu[i]))
				if(exp(mh.star.nu-mh.0.nu)>runif(1)){
					nu[i] <- nu.star[i]
					keep$nu[i] <- keep$nu[i]+1
					keep.tmp$nu[i] <- keep.tmp$nu[i]+1
				}
			}
		
			###
			###  Sample mu
			###
	 # browser()		
			
			mu.star[lc.idx,1] <- rnorm(T.tmp,mu.tmp[,1],tune$mu[i]) # proposals for mu
			mu.star[lc.idx,2] <- rnorm(T.tmp,mu.tmp[,2],tune$mu[i]) # proposals for mu

			idx.tmp <- lc.idx[which(!is.na(extract(S,mu.star[lc.idx,])))]
		
			s.tmp <- s[idx.tmp,]
			mu.tmp <- mu[idx.tmp,]
			
			mh.star.mu <- get.dmvt2(s.tmp,mu.star[idx.tmp,],sigma[i],a[i],rho[i],nu[i])
			mh.0.mu <- get.dmvt2(s.tmp,mu.tmp,sigma[i],a[i],rho[i],nu[i])
			keep.idx <- exp(mh.star.mu-mh.0.mu)>runif(length(idx.tmp))
			keep.idx <- idx.tmp[keep.idx]
			mu[keep.idx,] <- mu.star[keep.idx,]
			n.keep <- length(keep.idx)
			keep$mu[i] <- keep$mu[i]+n.keep
			keep.tmp$mu[i] <- keep.tmp$mu[i]+n.keep	
		}  # end loop though error classes (m)
			
			
		###
	 	###  Save samples
	 	###

		mu.save[,,k] <- mu
		sigma.save[k,] <- sigma
		a.save[k,] <- a
		rho.save[k,] <- rho
		nu.save[k,] <- nu
	}  # end loop through MCMC iterations (k)
	
	t.mcmc.end <- Sys.time()


	###
	###  Write Output 
	###
# browser()	
	keep <- lapply(keep, function(x) x/n.mcmc)
	keep$mu <- keep$mu/unlist(T.list)
	
	t.end <- Sys.time()
	cat(paste("\n\nEnd time:",t.end))
	cat(paste("\nTotal time elapsed:",round(difftime(t.end,t.start,units="hours"),2)),"hours")
	
	end <- list(mu=mu.save[,,k],sigma=sigma.save[k,],a=a.save[k,],rho=rho.save[k,],nu=nu.save[k,])
	list(mu=mu.save,sigma=sigma.save,a=a.save,rho=rho.save,nu=nu.save,
		keep=keep,end=end,tune=tune)
}