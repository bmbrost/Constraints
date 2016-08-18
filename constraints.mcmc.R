# Software Disclaimer:

# Although this software program has been used by the U.S. Geological Survey (USGS), no warranty, 
# expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and 
# functioning of the program and related program material nor shall the fact of distribution 
# constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.


constraints.MCMC <- function(s,S,lu,X,priors,start,tune,n.mcmc,n.cores=NULL,Cpp=TRUE){
	
	###
	### Brian M. Brost (09 MAR 2015)
	### 
	### See Brost et al. (Animal movement constraints improve resource selection inference in 
	### the presence of telemetry error), Appendix A for the model statement corresponding 
	### to this MCMC algorithm, and Appendix B for pseudocode pertaining to implementation.
	###
	
	t.start <- Sys.time()
	cat(paste("Start time:",t.start,"\n"))
		
	###
	### Libraries and Subroutines
	###
	
	library(mvtnorm) #For beta propsals
	library(doParallel) #For parallel processing
	library(foreach) #For parallel processing	
	library(data.table) #For fast database joins with the swim distance lookup table (lu)
	library(raster) #For spatial operations (e.g., using S as a mask for mu proposals)
	if(Cpp){ #For C++ routine to update mu
		library(Rcpp)
		library(inline)
		library(RcppArmadillo)
	}

	get.start <- function(mod,S){ #Get starting values from previous model output
		mod <- mod[[1]]
		k <- nrow(mod$sigma)
		mu.start <- mod$mu[k,]
		mu.start <- xyFromCell(S,mu.start)
		list(mu=mu.start,sigma=mod$sigma[k,],rho=mod$rho[k,],a=mod$a[k,],nu=mod$nu[k,],
			phi=mod$phi[k],beta=mod$beta[k,])
	}

	get.swimdist <- function(cell,lu,row.idx){ #Caclulate initial swim distances between consecutive locations
		N <- length(cell)
		cell.idx.tmp <- cbind(cell[-N],cell[-1])
		cell.idx.tmp <- t(apply(cell.idx.tmp,1,sort))
		d <- lu[J(cell.idx.tmp[,1],cell.idx.tmp[,2])]$dist
		d[is.na(d)] <- 0
		d <- c(NA,d)
		d[row.idx[,1]] <- NA
		d
	}
	
	get.dmvt2 <- function(x,y,sigma,a,r,nu){ #Calculate log density of mixture t-distribution (Eq. 1)
		d <- 2
		E.1 <- sigma^2*matrix(c(1,sqrt(a)*r,sqrt(a)*r,a),2,2,byrow=TRUE)
		E.2 <- sigma^2*matrix(c(1,sqrt(a)*-r,sqrt(a)*-r,a),2,2,byrow=TRUE)	
		c <- lgamma((nu+d)/2)-(lgamma(nu/2)+log(nu)+log(pi))+log(det(E.1)^-0.5)
		dif <- x-y
		kernel.1 <- (1+1/nu*(rowSums((dif%*%solve(E.1))*dif)))^-((nu+d)/2)
		kernel.2 <- (1+1/nu*(rowSums((dif%*%solve(E.2))*dif)))^-((nu+d)/2)
		log(0.5)+c+log((kernel.1+kernel.2))
	}

	
	###
	###  Create cluster for parallel processing
	###

	if(is.null(n.cores)) n.cores <- detectCores() 
	if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
	mcoptions <- list(preschedule=TRUE)
	cat(paste("\nUsing",n.cores,"cores for parallel processing."))

	###
	### Starting values: Appendix B, step 1
	###

	cat("\nGetting starting values....")
	if(length(start)==1){
		get.start(start,S)
		sigma <- get.start(start,S)$sigma
		a <- get.start(start,S)$a
		rho <- get.start(start,S)$rho
		mu <- get.start(start,S)$mu
		phi <- get.start(start,S)$phi
		beta <- get.start(start,S)$beta
		nu <- get.start(start,S)$nu
	}
	if(length(start)>1){
		sigma <- start$sigma
		a <- start$a
		rho <- start$rho
		mu <- start$mu
		nu <- start$nu
		phi <- start$phi
		beta <- start$beta
	}

	###
	###  Setup Variables 
	###
	
	N <- nrow(s) #Number of telemetry locations
	qX <- ncol(X) #Column dimension of design matrix
	m <- length(unique(s$lc)) #Number of error classes
	#Index of first and last location for each individual in data set
	row.idx <- t(sapply(unique(s$deploy.idx),function(x) range(which(s$deploy.idx==x))))
	#Index of locations in each error class
	lc.list <- sapply(sort(unique(s$lc)),function(x) which(s$lc==x),simplify=FALSE)
	
	# Location-level quantities
	lc <- as.numeric(s$lc) #Error class of each telemetry location s[t]
	cell <-  extract(S,mu,cellnumbers=TRUE)[,1] #Cell number of S in which mu[t] are located.
	cell.tilde <-  extract(S,mu) #Cell number of S.tilde in which mu[t] are located
	#Note: cell.tilde is used for matching mu[t] to the corresponding row in the design matrix X. 
	#When S and S.tilde are the same resolution, cell=cell.tilde. When inference pertaining to process
	#model parameters is at a coarser resolution than S, the support of mu, cell!=cell.tilde.
	delta.t <- s$delta.t #Time between consecutive locations: delta.t[t]=t[t]-t[t-1]
	delta.d <- get.swimdist(cell.tilde,lu,row.idx) #Distance between locations: delta.d[t]=d(mu[t-1],mu[t])
	
	s <- s@coords #Coordinates for observed telemetry locations s[t]
	lu <- copy(lu) #Lookup table of swim distances between all possible pairs of cells in S.tilde
	S.match <- sort(unique(values(S))) #Index for matching cell numbers to rows in X 
	Xbeta <- X%*%beta #Over all cells in S
	#Calculate integral for Appendix B, steps 3 and 4 (this quantitiy cancels in the mh ratio of step 5)
	int <- c(NA,foreach(x=2:N,.combine=c,.options.multicore=mcoptions) %dopar% 
		sum(exp(Xbeta-(lu[J(cell.tilde[x-1])]$dist/(delta.t[x]*phi))^1)))
	#Note: int[t]=integral{x(mu)*beta-d(mu,mu[t-1])/(delta.t[t]*phi)}dmu.
	#int is shifted such that int[t] is the integral wrt mu[t-1]; thus, int lines up 
	#with delta.t and delta.d, which simplifies calculations in the updates for phi, beta, and mu.

	int.star <- numeric(N) #Same as int, but for proposals (e.g., phi.star, beta.star, mu.star)
	delta.d.star <- matrix(0,N,2) #For mu.star: delta.d.star[t,]=c(d(mu[t-1],mu.star[t]),d(mu.star[t],mu[t+1]))
	X.match.star <- numeric(N) #Index for matching cell number of mu.star to rows in X
	Xbeta.mut.star <- numeric(N) #Resource selection component for mu.star: x(mu.star[t])*beta
			
	cat("\nCreating receptacles for output....")
	mu.save <- matrix(0,n.mcmc,N)
	sigma.save <- matrix(0,n.mcmc,m)
	a.save <- matrix(0,n.mcmc,m)
	rho.save <- matrix(0,n.mcmc,m)
	phi.save <- numeric(n.mcmc)
	beta.save <- matrix(0,n.mcmc,ncol(X))
	keep <- list(mu=rep(0,N),sigma=rep(0,m),a=rep(0,m),rho=rep(0,m),phi=0,beta=0,nu=rep(0,m))
	nu.save <- matrix(0,n.mcmc,m)


	###
	### Standardize parameters
	###

	cat("\nStandardizing variables....")
	s.scale <- max(apply(s,2,function(x) max(x)-min(x)))/6
	s.center <- apply(s,2,function(x) max(x)+min(x))/2
	d.scale <- max(lu$dist)/3
	d.center <- 0 # (max(lu$dist)-min(lu$dist))/2
	t.scale <- max(delta.t,na.rm=TRUE)/3

	s <- (s-matrix(s.center,nrow=N,ncol=2,byrow=TRUE))/s.scale
	mu <- (mu-matrix(s.center,nrow=N,ncol=2,byrow=TRUE))/s.scale
	sigma <- sigma/s.scale
	phi <- phi/d.scale*t.scale
	lu[,dist:=(dist-d.center)/d.scale]
	delta.d <- (delta.d-d.center)/d.scale
	delta.t <- delta.t/t.scale
	
	priors$sigma <- priors$sigma/s.scale
	priors$phi <- priors$phi/d.scale*t.scale
	tune$mu <- tune$mu/s.scale
	tune$sigma <- tune$sigma/s.scale
	tune$phi <- tune$phi/d.scale*t.scale
	tune$beta <- tune$beta*diag(qX) # tune$beta*solve(t(X)%*%X)


	###
	###  Begin MCMC Loop 
	###

	t.mu.update <- 0
	t.mcmc.start <- Sys.time()
	cat("\nEntering MCMC Loop....\n")

	for(k in 1:n.mcmc){ #Loop to iterate MCMC algorithm: Appendix B, step 7
		if(k%%100==0) cat(k,"");flush.console()
	
			###
			###  Sample components of Sigma (t=1:T): Appendix B, step 2
			###
	
			sigma.star <- rnorm(m,sigma,tune$sigma) #Proposals for sigma
			a.star <- rnorm(m,a,tune$a) #Proposals for a
			rho.star <- rnorm(m,rho,tune$rho) #Proposals for rho
			nu.star <- rnorm(m,nu,tune$nu) #Proposals for nu
					
			for(i in 1:m){ #Loop to iterate over error classes: Appendix B, step 2(f)
	
				idx <- lc.list[[i]] #Index of locations in error class i
	
				### Sample sigma: Appendix B, step 2(b)
	
				if(sigma.star[i]>priors$sigma[1] & sigma.star[i]<priors$sigma[2]){
					mh.star.sigma <- sum(get.dmvt2(s[idx,],mu[idx,],sigma.star[i],a[i],rho[i],nu[i]))
					mh.0.sigma <- sum(get.dmvt2(s[idx,],mu[idx,],sigma[i],a[i],rho[i],nu[i]))
					if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
						sigma[i] <- sigma.star[i]
						keep$sigma[i] <- keep$sigma[i]+1
					}
				}
	
				### Sample a: Appendix B, step 2(c)
	
				if(a.star[i]>priors$a[1] & a.star[i]<priors$a[2]){
					mh.star.a <- sum(get.dmvt2(s[idx,],mu[idx,],sigma[i],a.star[i],rho[i],nu[i]))
				  	mh.0.a <- sum(get.dmvt2(s[idx,],mu[idx,],sigma[i],a[i],rho[i],nu[i]))
					if(exp(mh.star.a-mh.0.a)>runif(1)){
						a[i] <- a.star[i]
						keep$a[i] <- keep$a[i]+1
					}
				}
	
				### Sample rho: Appendix B, step 2(d)
	
				if(rho.star[i]>priors$rho[1] & rho.star[i]<priors$rho[2]){
					mh.star.rho <- sum(get.dmvt2(s[idx,],mu[idx,],sigma[i],a[i],rho.star[i],nu[i]))
				  	mh.0.rho <- sum(get.dmvt2(s[idx,],mu[idx,],sigma[i],a[i],rho[i],nu[i]))
					if(exp(mh.star.rho-mh.0.rho)>runif(1)){
						rho[i] <- rho.star[i]
						keep$rho[i] <- keep$rho[i]+1
					}
				}
				
				### Sample nu: Appendix B, step 2(e)
	
				if(nu.star[i]>priors$nu[1] & nu.star[i]<priors$nu[2]){
					mh.star.nu <- sum(get.dmvt2(s[idx,],mu[idx,],sigma[i],a[i],rho[i],nu.star[i]))
			  		mh.0.nu <- sum(get.dmvt2(s[idx,],mu[idx,],sigma[i],a[i],rho[i],nu[i]))
					if(exp(mh.star.nu-mh.0.nu)>runif(1)){
						nu[i] <- nu.star[i]
						keep$nu[i] <- keep$nu[i]+1
					}
				}
			}
	
	
			###
			###  Sample phi: Appendix B, step 3
			###
	
			X.match <- match(cell.tilde,S.match) #Index for matching cell number of mu to rows in X
			phi.star <-  rnorm(1,phi,tune$phi) #Proposal for phi
			if(phi.star>priors$phi[1] & phi.star<priors$phi[2]){
				int.star <- c(NA,foreach(x=2:N,.combine=c,.options.multicore=mcoptions) %dopar% 
					sum(exp(Xbeta-(lu[J(cell.tilde[x-1])]$dist/(delta.t[x]*phi.star))^1)))
				mh.star.phi <- sum(X[X.match,]%*%beta-(delta.d/(delta.t*phi.star))^1-
					log(int.star),na.rm=TRUE)
				mh.0.phi <- sum(X[X.match,]%*%beta-(delta.d/(delta.t*phi))^1-
					log(int),na.rm=TRUE)
				if(exp(mh.star.phi-mh.0.phi) > runif(1)){
			    	phi <- phi.star   
					keep$phi <- keep$phi+1
					int <- int.star 
			 	}
			}
			
			
			###
			###  Sample beta: Appendix B, step 4
			###
	
			beta.star <- c(rmvnorm(1, beta, tune$beta)) #Proposal for beta
			Xbeta.star <- X%*%beta.star #Over all cells in S
			int.star <- c(NA,foreach(x=2:N,.combine=c,.options.multicore=mcoptions) %dopar% 
					sum(exp(Xbeta.star-(lu[J(cell.tilde[x-1])]$dist/(delta.t[x]*phi))^1)))
			mh.star.beta <- sum(X[X.match,]%*%beta.star-(delta.d/(delta.t*phi))^1-
					log(int.star),na.rm=TRUE) + dmvnorm(beta.star,priors$beta.mn,priors$beta.var,log=TRUE)
		  	mh.0.beta <- sum(X[X.match,]%*%beta-(delta.d/(delta.t*phi))^1-
		  		log(int),na.rm=TRUE) + dmvnorm(beta,priors$beta.mn,priors$beta.var,log=TRUE)
			if(exp(mh.star.beta-mh.0.beta) > runif(1)){
		    	beta <- beta.star   
	     		keep$beta <- keep$beta+1
				int <- int.star 
		  		Xbeta <- Xbeta.star
		  	}		
			
			
			###
			###  Sample mu: Appendix B, step 5
			###
	 
	 		t.mu.start <- Sys.time()
			mu.star <- matrix(rnorm(N*2,mu,tune$mu[lc]),N,2) #Proposals for mu
			cell.star <- rbind(c(NA,NA), extract(S,mu.star[-c(1,N),]*s.scale +
				matrix(s.center,N-2,2,byrow=TRUE),cellnumbers=TRUE),c(NA,NA)) #Cell numbers of mu.star
			idx <- which(!is.na(cell.star[,2]) & cell.star[,1]!=cell) #Index of proposals in S and different from k-1
			mu.star[idx,] <- (xyFromCell(S,cell.star[idx,1])-
				matrix(s.center,nrow=length(idx),ncol=2,byrow=TRUE))/s.scale #Standardize mu.star
			int.star[idx+1] <- foreach(x=idx,.combine=c,.options.multicore=mcoptions) %dopar% 	
				sum(exp(Xbeta-(lu[J(cell.star[x,2])]$dist/(delta.t[x+1]*phi))^1))
			tmp1 <- cbind(cell.tilde[idx-1],cell.star[idx,2]) #tmp1[t,]=c(cell(mu[t-1]), cell(mu.star[t]))
			tmp2 <- cbind(cell.tilde[idx+1],cell.star[idx,2]) #tmp2[t,]=c(cell(mu[t+1]), cell(mu.star[t])
			delta.d.star[idx,] <- cbind(lu[J(tmp1[,1],tmp1[,2])]$dist,
				lu[J(tmp2[,1],tmp2[,2])]$dist) #delta.d.star[t,]=c(d(mu[t-1],mu.star[t]), d(mu[t+1],mu.star[t]))
	
			X.match.star[idx] <- match(cell.star[idx,2],S.match) # Index for matching cell numbers to rows in X			
			Xbeta.mut.star[idx] <- X[X.match.star[idx],]%*%beta #Resource selection for mu.star: x(mu.star[t])*beta
			Xbeta.mut <- X[X.match,]%*%beta #Resource selection for mu: x(mu[t])*beta
			
			if(Cpp){ #Update mu with C++ routine
				mu.out <- updateMuCpp(idx,s,mu,mu.star,sigma[lc],a[lc],rho[lc],nu[lc],
					int,int.star,Xbeta.mut,Xbeta.mut.star,delta.d,delta.d.star,delta.t,phi)
				delta.d <- mu.out$delta.d #Update delta.d
				idx.keep <- mu.out$out==1 #Index of propsals (mu.star) to keep
			} else{ #Update mu with loop in R
				idx.keep <- logical(N)
				for(i in idx){
					mh.star.mu <- get.dmvt2(s[i,],mu.star[i,],sigma[lc[i]],a[lc[i]],rho[lc[i]],nu[lc[i]]) +
						Xbeta.mut.star[i] - (delta.d.star[i,1]/(delta.t[i]*phi))^1 - 
				   		(delta.d.star[i,2]/(delta.t[i+1]*phi))^1 - log(int.star[i+1])
				   	mh.0.mu <- get.dmvt2(s[i,],mu[i,],sigma[lc[i]],a[lc[i]],rho[lc[i]],nu[lc[i]]) +
				   		Xbeta.mut[i] - (delta.d[i]/(delta.t[i]*phi))^1 -
				   		(delta.d[i+1]/(delta.t[i+1]*phi))^1 - log(int[i+1])
					if(exp(mh.star.mu-mh.0.mu)>runif(1)){
						idx.keep[i] <- TRUE #Index of proposals (mu.star) t0 keep
						delta.d[c(i,i+1)] <- delta.d.star[i,] #Update delta.d
					}
				}
			}	
			idx <- which(idx.keep) #Index of proposals (mu.star) to keep
			mu[idx,] <- mu.star[idx,] #Update mu
			keep$mu[idx] <- keep$mu[idx]+1
			cell[idx] <- cell.star[idx,1] #Update cell number of S in which mu[t] is located
			cell.tilde[idx] <- cell.star[idx,2] #Update cell number of S.tilde which mu[t] is located
			int[idx+1] <- int.star[idx+1] #Update integral
	 	 			
			t.mu.update <- t.mu.update+difftime(Sys.time(),t.mu.start,"secs")
	
	
			###
		 	###  Save Samples: Appendix B, step 6
		 	###
	
			mu.save[k,] <- cell
			sigma.save[k,] <- sigma*s.scale
			a.save[k,] <- a
			rho.save[k,] <- rho
			phi.save[k] <- phi*d.scale/t.scale
			beta.save[k,] <- beta
			nu.save[k,] <- nu
	}
	t.mcmc.end <- Sys.time()
	registerDoSEQ()


	###
	###  Write Output 
	###
	
	keep <- lapply(keep, function(x) x/n.mcmc)
	keep$mu <- tapply(keep$mu,lc,mean)
	
	t.end <- Sys.time()
	cat(paste("\n\nEnd time:",t.end))
	cat(paste("\nTotal time elapsed:",round(difftime(t.end,t.start,units="hours"),2)),"hours")
	cat(paste("\nTime per MCMC iteration:",round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,
		2)),"seconds")
	cat(paste("\nTime per mu update:",round(t.mu.update/n.mcmc,2),"seconds"))
	system("say I am finally done")
	
	list(mu=mu.save,sigma=sigma.save,a=a.save,rho=rho.save,phi=phi.save,beta=beta.save,nu=nu.save,keep=keep)
}