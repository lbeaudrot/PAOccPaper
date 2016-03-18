#Test code

data.generator<-function(points,days,psi,p,phi,gamma,years) {
        require(unmarked)
        #first year of data
        y1<-matrix(NA,nr=points,nc=days)
        #generate the expected occupancies
        z<-rbinom(points,1,psi)
        #generate the observations
        for(i in 1:points)
                y1[i,]<-rbinom(days,1,z[i]*p)
        #subsequent years
        #three dimensional matrix to store the results
        yk<-array(NA,dim=c(points,days,years))
        yk[,,1]<-y1
        for(k in 2:years){
                #generate the deterministic part of the model
                occ<-apply(yk[,,k-1],1,max,na.rm=T)
                z<-rbinom(points,1,occ*phi+(1-occ)*gamma)
                #generate the observations
                for(i in 1:points)
                        yk[i,,k]<-rbinom(days,1,z[i]*p)
                
        }  
        # convert results to a two dimensional matrix for colext
        yk <- matrix(yk, points, days*years)
        #ny <- matrix(as.character(1:years),nrow(yk), years, byrow=TRUE)
        # store data in an unmarkedMultFrame object so it is ready for analysis
        yUMF <- unmarkedMultFrame(y = yk, numPrimary = years)
        yUMF
} 

set.seed(400)
data <- data.generator(points = 30,days = 30, psi = 0.5, p = 0.2, phi = 0.2, 
                       gamma = 0.1, years = 5)

#Look at the first year of data
summary(data)
model <- colext(psiformula = ~1,gammaformula = ~1,epsilonformula = ~1, 
                pformula = ~1, data = data, method = "BFGS",se = FALSE)
summary(model)
backTransform(model,type=c("psi"))
backTransform(model,type=c("col"))
backTransform(model,type=c("ext"))
backTransform(model,type=c("det"))

data3d <- array(data@y,dim = c(30,30,5))
obs <- apply(apply(data3d,c(1,3),function(z){max(z)}),2,sum)/30
mod <- as.numeric(smoothed(model)[,2])

plot(1:5, obs, type="p")
lines(1:5, mod)
