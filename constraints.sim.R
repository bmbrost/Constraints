rm(list=ls())
library(sp)
library(raster)
library(rgeos)
library(rstan)
# library(adehabitatHR)

setwd("/Users/Brost/Documents/research/publications/environmetrics/population_inference/sim/")
load("sim.Rdata")
# save.image("sim.Rdata")
ls()


#######################################################################
### Simulate data
#######################################################################

###
### Simulate environmental covariate for resource selection
###

# Spatial support of point process
S <- raster(xmn=-0.1,xmx=1.1,ymn=-0.1,ymx=1.1,res=0.02)

# Random resource covariate
d <- as.matrix(dist(xyFromCell(S,1:ncell(S))))		
phi <- 0.05
Sigma <- exp(-d/phi)  
Sigma <- t(chol(Sigma))
values(S) <- Sigma %*% rnorm(ncell(S))
plot(S)

# Rescale to finer resolution
S <- disaggregate(S,fact=2,method="bilinear")
S <- crop(S,raster(xmn=0,xmx=1,ymn=0,ymx=1))
plot(S)

# Center and scale covariates
values(S) <- scale(values(S))  # center and scale covariate

X <- cbind(1,values(S))


###
### Simulate animal locations
###

p <- nlayers(S)+1  # number of covariates
J <- 20  # number of individuals
X.list=vector("list",J)
y.list=vector("list",J)
s.list <- vector("list",J)
beta.list=vector("list",J)
tau.list=vector("list",J)
# hr <- cbind(runif(J,0.05,0.95),runif(J,0.05,0.95))  # homerange centroid
hr <- matrix(0,J,2)

mu.beta <- c(0,1)  # population mean of beta
s2.beta <- 0.25
n <- 1000  # number of proposed locations
for(j in 1:J){
# j <- 1
	# Simulate individual-level parameters
	beta.list[[j]] <- rnorm(p,mu.beta,sqrt(s2.beta))
	
	# Simulate 'homerange'
	hr[j,] <- xyFromCell(S,sample(1:ncell(S),1,prob=exp(X%*%beta.list[[j]])))
	# plot(S[[1]]);points(hr[j,1],hr[j,2])
	hr.tmp <- gBuffer(SpatialPoints(matrix(hr[j,],1,2)),width=runif(1,0.15,0.25))
	# plot(hr.tmp,add=TRUE)	

	# Simulate locations	
	s.tmp <- cbind(runif(n,xmin(S),xmax(S)),runif(n,ymin(S),ymax(S)))
	s.tmp <- gIntersection(SpatialPoints(s.tmp),hr.tmp)  # reject locations not in hr
	X.tmp <- cbind(1,raster::extract(S,s.tmp))
	lambda <- exp(X.tmp%*%beta.list[[j]])  # intensity of point process
	keep <- rbinom(nrow(s.tmp@coords),1,prob=lambda/max(lambda))  # keepers
	s.tmp <- s.tmp[keep==1,]  # used locations
	s.list[[j]] <- s.tmp
	# points(s.tmp)

	# Create UD
	# ud <- kernelUD(s.list[[1]])
	# ver <- getverticeshr(test,99.9)
	# plot(ver,add=TRUE)

	# Create design matrix
	X.tmp <- raster::extract(S,hr.tmp,cellnumbers=TRUE)[[1]]
	X.list[[j]] <- cbind(1,X.tmp[,2])

	# Create response variable
	cell.idx <- raster::extract(S,s.tmp,cellnumbers=TRUE)[,1]
	cell.tab <- table(cell.idx)
	y.tmp <- numeric(nrow(X.tmp))
	m <- match(as.numeric(names(cell.tab)),X.tmp[,1])
	y.tmp[m[which(!is.na(m))]] <- cell.tab[which(!is.na(m))]
	y.list[[j]] <- y.tmp
	# points(xyFromCell(S,X.tmp[which(y.tmp==1),1]),col=2,cex=0.5)

	tau.list[[j]]=rep(1,length(y.tmp))	
}
rm(cell.idx,cell.tab,hr.tmp,j,keep,lambda,m,s.tmp,X.tmp,y.tmp)

n <- unlist(lapply(s.list,function(x) nrow(x@coords)))
plot(S)
lapply(s.list,points)

# Reformat beta.list
beta.mat <- matrix(unlist(beta.list),p)


#######################################################################
### Plot simulated data
#######################################################################

ramp.gray <- colorRampPalette(c("white","gray30"))(2000)  # Color palette

s <- do.call(rbind,lapply(s.list,function(x) x@coords))
idx <- unlist(sapply(1:J,function(x) rep(x,n[x])))  # individual idx

# Create inset map
# locator()
inset <- c(0.53,0.39)
inset <- gBuffer(SpatialPoints(matrix(inset,1,2)),width=0.075)
inset <- as(extent(as.vector(t(bbox(inset)))),"SpatialPolygons")
plot(inset,add=TRUE)

pdf("sim.data.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(1,1,1,1)) # par(mar=c(0,0,0,0))  
image(S,col=ramp.gray,asp=TRUE,axes=FALSE,xlim=c(0,1),ylim=c(0,1))
points(s,pch=19,col=rgb(0,0,0,0.25),cex=0.5)
legend("topleft",legend=c("a.)"),bty="n",adj=c(1,0.5))
plot(inset,add=TRUE,border="white")
box()
image(S,col=ramp.gray,asp=TRUE,axes=FALSE,xlim=c(bbox(inset)[1,]),
	ylim=c(bbox(inset)[2,]))
points(s[which(idx==1),],pch=1,col=rgb(0,0,0,0.5))
points(s[which(idx==2),],pch=2,col=rgb(0,0,0,0.5))
points(s[which(idx==3),],pch=3,col=rgb(0,0,0,0.5))
points(s[which(idx==5),],pch=4,col=rgb(0,0,0,0.5))
points(s[which(idx==11),],pch=5,col=rgb(0,0,0,0.5))
points(s[which(idx==13),],pch=6,col=rgb(0,0,0,0.5))
points(s[which(idx==20),],pch=22,col=rgb(0,0,0,0.5))
# text(s,label=as.character(idx),col=rgb(0,0,0,0.5),cex=1)
legend("topleft",legend=c("b.)"),bty="n",adj=c(1,0.5))
box()
dev.off()


#######################################################################
### Fit hierarchical Poisson-normal Model 
#######################################################################

# Fit the independent Poisson model below to get reasonable tuning parameter for beta
source("code/pois.N.hier.mcmc.R")
a <- Sys.time()
mcmc.1.out=pois.N.hier.mcmc(y.list,X.list,20000,.32,FALSE,tau.list)
Sys.time()-a

layout(matrix(1:p,p,1))
matplot(t(mcmc.1.out$beta.save[1,,]),type="l",lty=1,col=rgb(0,0,0,alpha=.3))
lines(mcmc.1.out$mu.save[,1],col=2)
abline(h=mu.beta[1],col=8,lwd=2)
matplot(t(mcmc.1.out$beta.save[2,,]),type="l",lty=1,col=rgb(0,0,0,alpha=.3))
lines(mcmc.1.out$mu.save[,2],col=2)
abline(h=mu.beta[2],col=8,lwd=2)

plot(mcmc.1.out$beta.save[2,20,],type="l")

# Effective sample sizes
effectiveSize(mcmc(mcmc.1.out$mu.save[,2]))  # mu.beta
mean(sapply(1:J,function(x) effectiveSize(mcmc(mcmc.1.out$beta.save[2,x,]))))


#######################################################################
### Fit independent models and Lunn hierarchical model
#######################################################################

###
### Fit independent Poisson-normal models in parallel
###

# Use adaptive tuning to get tuning parameter for hierarchical model above
cl <- makeCluster(8)  # create cluster
registerDoParallel(cl)  # register cluster
mcoptions <- list(preschedule=TRUE)
source('code/pois.N.hier.mcmc.R')
a <- Sys.time()
mcmc.2.out=foreach(x=1:J,.options.multicore=mcoptions) %dopar%
	pois.N.hier.mcmc(y.list[x],X.list[x],20000,.32,TRUE,tau.list[x],adapt=TRUE)
Sys.time()-a
stopCluster(cl)  # kill parallel backend
hist(unlist(lapply(mcmc.2.out,function(x) x$keep)))
mean(unlist(lapply(mcmc.2.out,function(x) x$beta.tune)))

# # Examine output
idx <- 2
par(mfrow=c(2,1))
plot(mcmc.2.out[[idx]]$beta.save[1,,],type="l");abline(h=beta.mat[1,idx],col=2)
plot(mcmc.2.out[[idx]]$beta.save[2,,],type="l");abline(h=beta.mat[2,idx],col=2)

# Effective sample sizes
mean(unlist((lapply(mcmc.2.out,function(x) apply(x$beta.save[,,],1,function(y) effectiveSize(mcmc(y)))))))

###
### Fit independent Poisson-normal model in Stan
###

# # Prime Stan algorithm
# idx <- 1
# sim.stan.out=stan(file='code/rsf_model.stan',
	# data=list(N=length(y.list[[idx]]),y=y.list[[idx]],x1=X.list[[idx]][,2],p=2),
	# iter=100,chains=1)

# out.array=extract(sim.stan.out,permuted=FALSE,inc_warmup=TRUE)
# par(mfrow=c(p,1))
# plot(out.array[,,1],type="l");abline(h=beta.mat[1,idx],col=2)
# plot(out.array[,,2],type="l");abline(h=beta.mat[2,idx],col=2)

# mcmc.2.out <- vector("list",J)
# a <- Sys.time()
# for(i in 1:J){
	# mcmc.2.out[[i]]=stan(file='code/rsf_model.stan',
		# data=list(N=length(y.list[[i]]),y=y.list[[i]],x1=X.list[[i]][,2],p=2),
		# iter=20000,chains=1)
# }
# Sys.time()-a

# mcmc.2.out <- lapply(mcmc.2.out,function(x) extract(x,permuted=FALSE,inc_warmup=TRUE))

# idx <- 2
# plot(mcmc.2.out[[idx]][,,1],type="l");abline(h=beta.mat[1,idx],col=2)
# plot(mcmc.2.out[[idx]][,,2],type="l");abline(h=beta.mat[2,idx],col=2)
# mean(mcmc.2.out[[idx]][,,2])

# # Effective sample sizes
# mean(sapply(1:J,function(x) effectiveSize(mcmc(mcmc.2.out[[x]][,,1]))))


###
### Fit Lunn model
###

# Combine individual level parameters into same object for pois.N.Lunn.mcmc.R
beta.all <- array(0,c(p,J,20000))
for(i in 1:J){
	# beta.all[1,i,] <- mcmc.2.out[[i]][,,1]  # for STAN output
	# beta.all[2,i,] <- mcmc.2.out[[i]][,,2]  # for STAN output
	beta.all[1,i,] <- mcmc.2.out[[i]]$beta.save[1,,]  # for pois.N.hier.mcmc output
	beta.all[2,i,] <- mcmc.2.out[[i]]$beta.save[2,,]  # for pois.N.hier.mcmc output
}

# Check output; croscheck with output from mcmc.2.out to confirm correct transformation
idx <- 2
plot(beta.all[1,idx,],type="l");abline(h=beta.mat[1,idx],col=2)
plot(beta.all[2,idx,],type="l");abline(h=beta.mat[2,idx],col=2)

# Create model object for pois.N.Lunn.mcmc.R
out <- list(beta.save=beta.all,s2.init=100,mu.init=rep(0,p))

source("code/pois.N.Lunn.mcmc.R")
a <- Sys.time()
mcmc.3.out=pois.N.Lunn.mcmc(out)
Sys.time()-a

layout(matrix(1:p,p,1))
matplot(t(mcmc.3.out$beta.save[1,,]),type="l",lty=1,col=rgb(0,0,0,alpha=.3))
lines(mcmc.3.out$mu.save[,1],col=2)
abline(h=mu.beta[1],col=8,lwd=2)
matplot(t(mcmc.3.out$beta.save[2,,]),type="l",lty=1,col=rgb(0,0,0,alpha=.3))
lines(mcmc.3.out$mu.save[,2],col=2)
abline(h=mu.beta[2],col=8,lwd=2)

# Effective sample sizes
effectiveSize(mcmc(mcmc.3.out$mu.save[,2]))  # mu.beta
mean(apply(mcmc.3.out$beta.save[2,,],1,function(x) effectiveSize(mcmc(x))))


#######################################################################
### Plot parameter estimates
#######################################################################

get.sum <-function(x){
	c(mean(x),quantile(x,c(0.025,0.25,0.75,0.975)))
}

# Summarize results
mu.beta.1 <- get.sum(mcmc.1.out$mu.save[,2])  # mu.beta from hierarchical model
beta.1 <- t(apply(mcmc.1.out$beta.save[2,,],1,get.sum))  # beta from hierarchical model
mu.beta.3 <- get.sum(mcmc.3.out$mu.save[,2])  # mu.beta from Lunn model
beta.3  <- t(apply(mcmc.3.out$beta.save[2,,],1,get.sum))  # beta from Lunn model

w <- 0.175
lwd.95 <- 0.75
lwd.50 <- 2.25
pt.cex <- 0.65
col.lunn <- "gray55"

pdf("sim.results.pdf",width=6,height=9)
par(mar=c(2.5,2.1,1.1,1.1))
plot(1,ylim=c(1-0.5,J+1+0.5),xlim=range(c(0,c(rbind(beta.1,beta.3))))+c(-0.1,0.1),
	pch="",axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
polygon(c(-10,10,10,-10,-10),c(J-0.5,J-0.5,J+1.5,J+1.5,J-0.5)+1,col="gray95")
points(mu.beta.1[1],J+1+w,pch=19,cex=pt.cex)
lines(mu.beta.1[c(3,4)],rep(J+1+w,2),lwd=lwd.50)
lines(mu.beta.1[c(2,5)],rep(J+1+w,2),lwd=lwd.95)
points(mu.beta.3[1],J+1-w,pch=19,cex=pt.cex,,col=col.lunn)
lines(mu.beta.3[c(3,4)],rep(J+1-w,2),lwd=lwd.50,col=col.lunn)
lines(mu.beta.3[c(2,5)],rep(J+1-w,2),lwd=lwd.95,col=col.lunn)
for(i in 1:J){
	points(beta.1[i,1],i+w,pch=19,cex=pt.cex)
	lines(beta.1[i,c(3,4)],rep(i+w,2),lwd=lwd.50)
	lines(beta.1[i,c(2,5)],rep(i+w,2),lwd=lwd.95)
	points(beta.3[i,1],i-w,pch=19,cex=pt.cex,col=col.lunn)
	lines(beta.3[i,c(3,4)],rep(i-w,2),lwd=lwd.50,col=col.lunn)
	lines(beta.3[i,c(2,5)],rep(i-w,2),lwd=lwd.95,col=col.lunn)
}
box()
abline(v=0,lty=2)
axis(1,at=seq(-10,10,1))
axis(1,at=seq(-10,10,0.5),labels=FALSE)
# axis(2,at=1:(J+1),labels=FALSE)
sapply(1:J,function(x) 
	axis(2,at=J+1-x,tick=FALSE,label=substitute(beta[i],list(i=x)),las=1,pos=-0.025)) 
axis(2,at=J+1,tick=FALSE,label=expression(mu[beta]),las=1,pos=-0.025)
dev.off()

mean(mcmc.1.out$mu.save[,2])
mean(mcmc.3.out$mu.save[,2])

abline(v=mean(mcmc.1.out$beta.save[2,19,]))
abline(v=mean(mcmc.3.out$beta.save[2,19,]))
plot(mcmc.1.out$beta.save[2,1,],type="l")
plot(mcmc.3.out$beta.save[2,1,],type="l")