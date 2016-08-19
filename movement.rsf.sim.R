###
### Simulate true and observed locations in seal domain with temporal dependence and resource selection
###	Simulated data corresponds to model movement.rsf.MCMC
###

rm(list=ls())
setwd("~/Documents/research/harbor_seals/projects/constraints")

#############################################################################################
### Libraries and functions
#############################################################################################

library(sp)
library(rgdal)
library(maptools)
library(raster)
library(gstat)
library(rgeos)
library(adehabitatHR)
library(splancs)
library(mvtnorm)
library(MCMCpack)
library(MASS)
library(lattice)
library(gdistance)
library(parallel)
library(doParallel)
library(foreach)	
library(FNN)
library(data.table)

library(Rcpp)
library(inline)
library(RcppArmadillo)
Rcpp::sourceCpp("update.mu.cpp") #C++ routine for mu update (Appendix B, step 5)

get.swimdist.lu <- function(S){ #Create lookup table of swim distances between all pairwise cells in S
	library(data.table)
	library(raster)
	tran <- transition(S,function(x) 1,directions=4)
	idx <- which(!is.na(values(S)))
	lu.tmp <- c(costDistance(tran,xyFromCell(S,idx)))*xres(S)
	lu <- expand.grid(V1=idx,V2=idx)[,2:1]
	same.idx <- which(lu$V1==lu$V2)
	dup.idx <- unlist(sapply(0:(length(idx)-2),function(x) (x+1)*(length(idx)+1)-seq(0,x,1),simplify=TRUE))
	order.idx <- order(lu[dup.idx,2])
	lu$dist <- 0
	lu$dist[-c(dup.idx,same.idx)] <- lu.tmp
	lu$dist[dup.idx][order.idx] <- lu.tmp
	lu <- data.table(lu,key=c("V1","V2"))
	lu
}

get.mu <- function(T=1000,mu1,S,lu,X,beta,phi,lc=3,a,s){ #Simulate true locations (mu[t])
	library(raster)
	library(data.table)
	
	# d.scale <- max(lu$dist)/3
	# lu[,dist:=dist/d.scale]
	# X[,2] <- X[,2]/d.scale
	# beta <- beta*d.scale

	mu <- matrix(NA,T,5);colnames(mu)<-c("x","y","dtime","deploy.idx","lc")
	mu <- as.data.frame(mu)
	mu$dtime <- 0
	mu[1,1:2] <- mu1@coords
	mu$deploy.idx <- 1

	for(j in 2:T){
		mu$dtime[j] <- mu$dtime[j-1]+rgamma(1,shape=a,scale=s)
		dt <- mu$dtime[j]-mu$dtime[j-1]
		cell.idx <- extract(S,mu[j-1,1:2])
		lu.tmp <- as.matrix(lu[J(cell.idx)])
  		prob <- exp(X%*%beta-((lu.tmp[,3])/(dt*phi))^1)
		prob <- prob[match(values(S),lu.tmp[,1],nomatch=0)]
		cell.idx <- sample(which(!is.na(values(S))),1,prob=prob)
		mu[j,1:2] <- xyFromCell(S,cell.idx)
	}

	mu$dtime <- as.POSIXct(mu$dtime*60^2,origin="1960-01-01") # Time
	coordinates(mu) <- ~x+y # Convert to SpatialPoints
	mu$lc <- sample(1:lc,length(mu),replace=TRUE)
	mu$delta.t <- as.numeric(unlist(tapply(mu$dtime,mu$deploy.idx,function(x) 
		c(NA,difftime(x[-1],x[-length(x)],units="hours"))))) #Units are hours
	mu$delta.t <- ifelse(mu$delta.t==0,0.1,mu$delta.t)
	mu	
}

get.s <- function(mu,sigma,a,rho,nu){ #Simulate observed locations (s[t])
	library(mvtnorm)

	s <- mu
	T <- nrow(s)
	lc <- sort(unique(s$lc))
	Sigma <- lapply(lc, function(x) sigma[x]^2*matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))
	tdist.idx <- sample(1:2,T,replace=TRUE) # Index for mixture t-distribution
	K <- matrix(c(-1,0,0,1),2)
	for(i in 1:T){
		lc <- s$lc[i]
		if(tdist.idx[i]==1){	
			s@coords[i,] <- mu@coords[i,]+rmvt(1,sigma=Sigma[[lc]],df=nu[lc])
		}
		if(tdist.idx[i]==2){
			s@coords[i,] <- mu@coords[i,]+rmvt(1,sigma=K%*%Sigma[[lc]]%*%t(K),df=nu[lc])
		}
	}
	s
}


#############################################################################################
### Prepare spatial data for analysis
#############################################################################################

###
### Define S (support of mu)
###

proj.aea <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 		+ellps=GRS80 +towgs84=0,0,0"
ak.simple <- readOGR(dsn="/Users/brost/Documents/research/harbor_seals/geodata/base/physical_features",
	layer="ak_simple_poly_R")
ak <- readOGR(dsn="/Users/brost/Documents/research/harbor_seals/geodata/base/physical_features", 
	layer="ak_63360_poly_R")


# Define S according to some bounding box
# S.clip <- matrix(c(55000,110000,110000,55000,55000,805000,805000,840000,840000,805000),5,2)
# S.clip <- SpatialPolygons(list(Polygons(list(Polygon(S.clip)),"1")),proj4string=CRS(proj.aea))
# S.poly <- gIntersection(ak.simple,S.clip)

# Define S according to haul out location or centroid
haulout <- matrix(c(61156,830265),1) #PV95KOD12
haulout <- SpatialPoints(haulout,proj4string=CRS(proj.aea))

# Create polygon of S
S.buffer <- gBuffer(haulout,width=65000)
S.poly <- gIntersection(ak,S.buffer,byid=TRUE)
plot(S.poly,col="grey")

# Create raster of S
S <- raster(S.poly,resolution=500)
S <- rasterize(S.poly,S,field=1)
S <- reclassify(S,matrix(c(1,NA,NA,1),2,2,byrow=TRUE))
plot(S)

# Confirm haul-out site is in S
extract(S,haulout)
idx <- which(!is.na(values(S)))
neigh <- get.knnx(xyFromCell(S,idx),haulout@coords,k=1)$nn.index
haulout <- SpatialPoints(xyFromCell(S,idx[neigh]),proj4string=CRS(proj.aea))
extract(S,haulout)

# Remove isolated water cells enforcing 4-neighbor rule
# Note: make sure haulout is inside S
tran <- transition(S,function(x) 1,directions=4)
idx <- which(values(S)==1)
cd <- costDistance(tran, haulout, xyFromCell(S,idx))*xres(S)
values(S)[idx[which(cd==Inf)]] <- NA

# Limit extent of S to some distance d from haul-out location
d <- 60000
values(S)[idx[which(cd>d)]] <- NA
plot(S)

###
### Aggregate S to coarser resolution for network distance calculations and inference on Beta
###

S.tilde <- S
agg.factor <- 1
# agg.factor <- 2
# S.tilde <- aggregate(S,fact=agg.factor,fun=min,expand=TRUE)

### Asign cells in S corresponding cell number in S.tilde
idx <- which(!is.na(values(S)))
values(S)[idx] <- extract(S.tilde,xyFromCell(S,idx),cellnumbers=TRUE)[,1]
plot(S)

plot(S.poly)

### Create data.table containing all possible distances over network defined by S.tilde
lu <- get.swimdist.lu(S.tilde)


#############################################################################################
### Create resource selection covariates and design matrix
#############################################################################################

###
### Distance to haul-out
###

d.haulout <- Which(S,na.rm=FALSE)
tran <- transition(d.haulout,function(x) 1,directions=4)
idx <- which(values(d.haulout)==1)
cd <- costDistance(tran, haulout, xyFromCell(d.haulout,idx))*xres(d.haulout)
d.haulout[idx] <- cd
d.haulout <- aggregate(d.haulout,fact=agg.factor,fun=mean,na.rm=TRUE,expand=TRUE) 
plot(d.haulout)


###
### Distance to shore
###

# d.shore <- raster(S.poly,resolution=250)
# d.shore <- rasterize(S.poly,d.shore,field=1)
# d.shore <- distance(d.shore)
# d.shore <- crop(d.shore,S)
# idx <- which(is.na(values(S)))
# d.shore[idx] <- NA
# plot(S)
# plot(d.shore,add=TRUE)


###
### Bathymetry
###

# bath <- readOGR(dsn="/Users/brost/Documents/research/harbor_seals/projects/constraints/bathymetry",
	# layer="bath")

# Remove point sounding not in water defined by raster
# idx <- !is.na(extract(d.haulout,bath))
# bath <- bath[idx,]

# plot(bath,pch="")
# plot(d.haulout,add=TRUE)
# plot(bath,add=TRUE,pch=19,cex=0.1)

# # Identify cells along shoreline
# shore <- raster(S.poly,resolution=xres(S))
# shore <- rasterize(S.poly,shore,field=1)
# shore <- distance(shore)
# shore <- Which(shore==xres(S))
# plot(shore)

# # Create point data set for shoreline
# idx <- which(values(shore)==1)
# shore.pts <- xyFromCell(shore,idx)
# shore.pts <- data.frame(shore.pts,survy_d="shoreline",depth=0,qlty_cd=NA,active=NA)
# shore.pts <- SpatialPointsDataFrame(shore.pts[,1:2], shore.pts[,3:6], proj4string=CRS(proj.aea))
# plot(shore.pts)

# # Combine bathymetry data with shoreline
# bath <- rbind(bath,shore.pts)
# plot(bath,add=TRUE,col=2)

# # Use IDW to interpolate depths to centroids of raster cells
# idx <- which(values(d.haulout)>=0)
# xy.idx <- xyFromCell(d.haulout,idx)
# xy.idx <- SpatialPoints(xy.idx,proj4string=CRS(proj.aea))
# bath.idw <- idw(depth~1,locations=bath,xy.idx,idp=2,nmax=10)

# # Create raster with interpolated depths
# bath.grd <- d.haulout
# bath.grd[idx] <- bath.idw$var1.pred
# plot(bath.grd)
# writeRaster(bath.grd,"/Users/brost/Documents/research/harbor_seals/projects/constraints/bathymetry/bathymetry.grd",overwrite=TRUE)

bath <- raster("/Users/brost/Documents/research/harbor_seals/projects/constraints/bathymetry/bathymetry.grd")
plot(bath,col=colorRampPalette(c("lightblue","mediumblue"))(200))

###
### Create design matrix and covariate for simulation
###

idx <- which(!is.na(values(d.haulout)))
X <- cbind(1,values(d.haulout)[idx],values(bath)[idx])
plot(X[,2:3],pch=19,cex=0.1)
cor(X[,2:3])

# Center and scale design matrix
X.center <- c(0,apply(X[,-1],2,mean))
X.sd <- c(1,apply(X[,-1],2,sd))
X.scale <- t(apply(X,1,function(x) (x-X.center)/X.sd))


#####
##### Remove unneseccary objects, save workspace
#####

rm(idx,cd,S.buffer,neigh)
ls()


#############################################################################################
### Simulate true and observed locations
#############################################################################################

# Set parameters for observation model
sigma <- c(2291,2727,13252)
nu <- c(16.74,1.60,1.00)
a <- c(0.70,0.50,0.75)
rho <- c(0.85,0.16,0.30)

# Set parameters for process model
phi <- 424 #m/hr
beta <- c(0,-0.0001511532, -0.0270482365) #Not standardized

# Explore gamma distribution for simulating delta.t, the time between consecutive locations
# gamma.mm <-function(mu,sd){
	# a <- (mu^2)/(sd^2)
	# s <- (sd^2)/mu
	# list(shape=a,scale=s)
# }

# gamma.mm(4.3,6.1)
# hist(rgamma(10000,shape=1.1,scale=4.5),breaks=100)

# Simulate true locations with resource selection and temporal dependence
mu <- get.mu(T=300,haulout,S,lu,X,beta,phi,lc=3,a=1.1,s=4.5) #True locations

# Simulate observed locations
s <- get.s(mu,sigma,a,rho,nu) #Observed locations

plot(mu@coords,type="l",cex=log(mu$delta.t))
plot(S,add=TRUE)
points(mu[,1:2],type="b",cex=log(mu$delta.t))

plot(s[,1:2],pch="")
plot(S,add=TRUE)
points(s[,1:2],col=(s$lc+1))
points(mu[,1:2],col=1,pch=19,cex=0.5)


#############################################################################################
### Save/load workspace
#############################################################################################

# rm(list=ls())
setwd("~/Documents/research/harbor_seals/projects/constraints")
# save.image('~/Dropbox/constraints/figures/Figure_1.Rdata')
# save.image('~/Dropbox/constraints/supplement/workspace.Rdata')
# save.image("sim/sim.workspace.Rdata")
# load("sim/sim.workspace.Rdata")


#############################################################################################
### Fit model
#############################################################################################

#Prior distribution parameters
priors <- list(sigma=c(0,20000),phi=c(0,750),a=c(0,1),rho=c(0,1),nu=c(0,30),beta.mn=rep(0,3),beta.var=diag(3)*10^2)

#Set tuning parameters
tune <- list(mu=c(950,1200,1250),sigma=c(400,750,2750),rho=c(0.1,0.39,0.44),a=c(0.205,0.25,0.275),
	nu=c(15.0,0.59,0.4),phi=90,beta=0.045)

#Set starting values
start <- list(mu=mu@coords,sigma=sigma,nu=nu,a=a,rho=rho,phi=phi,beta=beta*X.sd)

source("movement.rsf.MCMC.R")
# start <- list(out1)
out2 <- movement.rsf.MCMC(s,S,lu,X.scale,priors,start,tune,n.mcmc=1000,n.cores=7,Cpp=TRUE)	

mod <- out2
mod$keep
matplot(mod$sigma,type="l",lty=1);abline(h=sigma,col=1:5,lty=3)
matplot(mod$a,type="l",lty=1);abline(h=a,col=1:5,lty=3)
matplot(mod$rho,type="l",lty=1);abline(h=rho,col=1:5,lty=3)
plot(mod$phi,type="l",lty=1);abline(h=phi,lty=3)
matplot(mod$beta[,-1],type="l",lty=1);abline(h=beta*X.sd,col=0:5,lty=3)
matplot(mod$nu,type="l",lty=1);abline(h=nu,col=1:5,lty=3)

matplot(rbind(out1$beta[,-1],out2$beta[,-1]),type="l",lty=1);abline(h=beta*X.sd,col=0:5,lty=3)
abline(v=5000)

apply(mod$sigma[-(1:100),],2,mean)
apply(mod$beta,2,mean)
apply(mod$nu,2,mean)
apply(mod$a,2,mean)
apply(mod$rho,2,mean)
mean(mod$phi)

apply(mod$beta,2,quantile,c(0.025,0.975))
apply(mod$sigma,2,quantile,c(0.025,0.975))

str(mod$mu)
plot(xyFromCell(S,mod$mu[,10])+rnorm(2000,0,30),type="p",col=colorRampPalette(c("green","red"))(1000),pch=19,cex=0.5)

points(xyFromCell(S,test1[10]),pch=19)


get.mode <- function(x){ # Find mode
	ux <- unique(x)
	ux[which.max(tabulate(match(x,ux)))]
}


test1 <- apply(mod$mu[,],2,get.mode)
test2 <- apply(out10$mu,2,get.mode)

test <- sapply(1:200, function(x) sqrt((xyFromCell(S,test1[x])[1]-xyFromCell(S,test2[x])[1])^2+
	(xyFromCell(S,test1[x])[2]-xyFromCell(S,test2[x])[2])^2))
summary(test)

test1 <- apply(out2$mu,2,get.mode)
test2 <- apply(out200$mu,2,get.mode)
table(test1==test2)







tab <- data.frame(table(c(mod$mu)))
post.mu <- S
idx <- which(!is.na(values(S)))
values(post.mu)[idx] <- 0
values(post.mu)[idx][match(tab[,1],idx)] <- tab[,2]
plot(post.mu)

mu.mode <- apply(mod$mu,2,get.mode)
mu.hat <- xyFromCell(S,mu.mode)
plot(mu.hat,pch="")
plot(post.mu,add=TRUE)
points(mu.hat,pch=19,cex=0.25,col=rgb(1,0,0))

plot(mu.hat,pch=19,cex=0.25,col=rgb(1,0,0))

dev.new()
plot(post.mu)
zoom(post.mu)

dev.new()
idx <- sample(1:nrow(mod$mu),1000)
plot(xyFromCell(S,c(mod$mu[idx,])),cex=0.05,col=rgb(1,0,0,0.01))
points(xyFromCell(S.prime,1:ncell(S.prime)),cex=0.1,pch=19)
points(mu.hat,pch=19,cex=0.1,col=1)
