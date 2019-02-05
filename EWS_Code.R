

#Gaurav Baruah, PhD student, University of Zurich, www.popecol.org; 
#https://gauravkbaruah.github.io/ ; (gaurav.baruah@ieu.uzh.ch)

rm(list=ls())
set.seed(12)

#environment function
env<-function(t){
  e.mean<-matrix(NA,nrow=5000,ncol=2)
  e.mean[1,1]<-0
  e.mean[1,2]<-0
  g<-0.0
  
  if (t <500) {
    
    e.mean[t,1]<-0
    e.mean[t,2]<-0
  }
  else if (t>=500) {
    e.mean[t,1]<-5 -5*(t/500)  # environmental cue
    e.mean[t,2]<- 20 -20*(t/500)     # optimum environment
  }          
 
  return(e.mean[t,])
}


#function for simulating the dynamics of the model
#dynamics() <- function that simulates the stochastic model
#tmax<- total time
#ind<- no. of individuals
#omega<- width of fitness function
#var.size<- variation in the mean.breeding value
#rho <- environmental predictability p
#var.U<- variability in cue
#var.theta <- variability in optimum phenotype
#r0 <- net reproductive rate.
#mean.breed<- mean breeding value

dynamics<-function(tmax,ind,omega,var.size,rho,var.U, var.Theta,r0,mean.breed,beta){
  breed<-matrix(NA,nrow=ind,ncol=tmax)
  W<-matrix(NA,nrow=ind,ncol=tmax)
  s<-matrix(NA,nrow=ind,ncol=tmax)
  c<-matrix(NA,nrow=tmax,ncol=2)
  mean.breed<-numeric()
  theta<-numeric()
  N<-numeric()
  pop.breed<-numeric()
  mean.fitness<-numeric()
  cue<-numeric()
  var<-numeric()
  phenotype<-numeric()
  rand.s<-numeric()
  
  #initialization
  var[1] <- var.size
  breed[,1]<-rnorm(ind,mean=0,sd=sqrt(var[1]))
  theta[1]<-0
  N[1]<-70
  mean.breed[1]<-phenotype[1] <-cue[1]<-rand.s[1]<-c[1,1]<-0
  pop.breed[1]<-mean(breed[,1])
  s[,1]<-breed[,1] 
  mean.fitness[1]<-1
  c[1,2]<-theta[1]

  
  for (t in 2:tmax) {
    
    S = matrix(c(var.U, rho*sqrt(var.U*var.Theta), rho*sqrt(var.U*var.Theta), var.Theta),2)  #S is a covariance matrix of the environment and the cue of the environment
    c[t,] <- mvrnorm(1,env(t),S)  #multivariate random distribution
    cue[t]<- c[t,1]# env(t)[1] +rnorm(1,0,sqrt(var.U )) # c[t,1]
    theta[t] <-  c[t,2] #env(t)[2] + rnorm(1,0,sqrt(var.Theta))#c[t,2]
    
    
    
    for (j in 1:ind){
      
      breed[j,t]<- rnorm(1,mean=mean.breed[t-1],sd=sqrt(var.size))  #population of breeding values
      pop.breed[t]<-mean(breed[,t]) 
      rand.s[t]<-rnorm(1,0,0.01)
      s[j,t] <- breed[j,t]+beta*cue[t]+rand.s[t]   #phenotypic value for individual j
      W[j,t] <-exp(-(s[j,t]-theta[t])^2/(2*omega)) #gaussian fitness function
      
    }
    
    
    mean.fitness[t] = mean(W[,t]) #mean fitness 
    var[t] <- var.size 
   
    #(var.size/(omega+var(s[,t-1])))*(theta[t]-pop.breed[t]-beta*cue[t]+rand.s[t]) 
    mean.breed[t]<-mean.breed[t-1]+ var.size*sqrt(omega)*(1/(var.size+omega))^1.5*(theta[t]-pop.breed[t]-beta*cue[t]+rand.s[t])*exp(-(pop.breed[t]+beta*cue[t]+rand.s[t]-theta[t])^2/(2*(var.size+omega)))
    phenotype[t]<- mean.breed[t] +beta*cue[t]+rand.s[t]  # dynamics of the phenotype
    N[t]<-N[t-1]*(r0^(1-N[t-1]/73))*mean.fitness[t]  # discrete population dynamics
    
  }
  
  list(N=N,mean.breed=mean.breed,phenotype=phenotype,theta=theta,cue=cue,one=mean.fitness,omega=omega,var=var)
  
}



### FOR parallelising the R code ####

mc.replica = function(iter, mean.breed, ...){
  set.seed(as.numeric(Sys.time())-Sys.getpid()) 
  init = mean.breed 
  replicate = dynamics(mean.breed=init, ...)
  replicate$start = init
  replicate$iter = iter
  return(replicate)
}
require(MASS)

#######################STARTING PARAMETERS##############

start.mean.breed=0
reps<-2 #no. of replicate simulation
EWS.trait.genvar.0.5<-lapply(1:reps, mc.replica,tmax=800,ind=200,omega=20,beta=0.2,var.size=0.5,rho=1,var.U=1.25^2, var.Theta=1.25^2,r0=1.2,mean.breed=start.mean.breed)

par(mfrow=c(2,2))
plot(EWS.trait.genvar.0.5[[1]]$N[500:530], typ='l')
s

