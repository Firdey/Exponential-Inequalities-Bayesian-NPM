library(fda)
library(mvtnorm)
library(MHadaptive)
sin_basis<-function(x,n){
  return(sqrt(2)*sin(2*pi*n*x));  
}
#Sin Basis Function

fourier_eval=function(theta,x){
  #Vector of theta k size and x n size vector.
  k=length(theta);
  eval = rep(0,length(x));
  eval = theta[1];
  #First coeff of basis is constant
  for (j in 2:k){
    if (j%%2==0) {
      eval=eval +theta[j]*sin_basis(x,j/2);
      #If even then a sin basis, sin basis is always firs  
    }
    else {
      eval=eval +theta[j]*cos_basis(x,(j-1)/2);
      #If odd then a cos basis, the index is shifted accordingly.
    }
  }
  return(eval)
}
cos_basis<-function(x,n){
  return(sqrt(2)*cos(2*pi*n*x));   
}
#Cos Basis Function

#f_0truth<-function(x){
# f_0=2+1.2*sin_basis(x,1)+1*cos_basis(x,1)+0.8*sin_basis(x,2)+0.7*cos_basis(x,2);
#  return(f_0)
#}
#Truth function we aim to be near to.

SimulationBM<-function(k,alpha,tao,n,sigma,Model){
  #Simulation of the quantity of interest with inputs:
  # k the number of components (would like to be greater than number of components of the f_0)
  # alpha, tao and sigma control
  # n the number of data
  
  x<-runif(n,0,1)
  # Generate the covariates (n of these)
  
  #Ignore
  #thetaPri=rep(0,k);
  #for (j in 1:k) {
  #theta[j]=rnorm(1,0,sqrt((tao^2)*j^(-2*alpha-1)));
  #}
  
  epsilon=rnorm(n,0,sigma);
  #Generate the errors with sigma 
  
  fourier_eval=function(theta,x){
    #Vector of theta k size and x n size vector.
    k=length(theta);
    eval = rep(0,length(x));
    
    for (j in 1:k){
      if (j==1) {
        eval = eval + theta[j];
        #First coeff of basis is constant
      }
      else if (j%%2==0) {
        eval=eval +theta[j]*sin_basis(x,j/2);
        #If even then a sin basis, sin basis is always first
        
      }
      else {
        eval=eval +theta[j]*cos_basis(x,(j-1)/2);
        #If odd then a cos basis, the index is shifted accordingly.
      }
    }
    return(eval)
  }
  
  if (Model==0){
    f_0truth<-function(x){
      f_0=2+1.2*sin_basis(x,1)+1*cos_basis(x,1)+0.8*sin_basis(x,2)+0.5*cos_basis(x,2);
      return(f_0)
    }
  }
  else if (Model==1){
    f_0truth<-function(x){
      res <- fourier_eval(coeff,x)
      return(res); 
    }
  }
  else if (Model==2){
    f_0truth<-function(x){
      res <- fourier_eval(coeff,x)
      return(res);
    }
  }
  else {
    f_0truth<-function(x){
      res <- fourier_eval(coeff,x)
      return(res);
    }
  }
  f_0truth_at_x<-f_0truth(x);
  y_data<-f_0truth_at_x + epsilon;
  #Create our observed noise data.
  
  PriorMult=function(theta,tao,alpha,k){
    PriorMult=0;
    for (l in 1:k){
      PriorMult=PriorMult+log(dnorm(theta[l],0,sqrt((tao^2)*l^(-2*alpha-1))));
    }
    return((PriorMult))
  }
  #The prior density multiplied as a function of theta implicitly.
  
  #Evaluates the value of f given theta and x. Note the vectorised format, output a vector of size n.
  Likelihood=function(y_data,x,epsilon,sigma,n,f_0truth_at_x,theta){
    likelihood=0;
    for (s in 1:n){
      likelihood=likelihood+log(dnorm(y_data[s],fourier_eval(theta,x)[s],sigma));
    }
    return((likelihood))
  }
  #Evaluates the likelihood, implicitly a function of theta.
  
  PosteriorKernel=function(theta,y_data,x,epsilon,sigma,n,f_0truth_at_x,tao,alpha,k){
    S=(Likelihood(y_data,x,epsilon,sigma,n,f_0truth_at_x,theta))+(PriorMult(theta,tao,alpha,k));
    #Obtain posterior Kernel without NC possibly.
    return(S)
  }
  
  Answer=(Metro_Hastings(function(theta) {PosteriorKernel(theta,y_data,x,epsilon,sigma,n,f_0truth_at_x,tao,alpha,k)},c(1.7,1,1,0.8,0.6),iterations = 500*200+4000,burn_in = 4000))
  return(list(mcmc_thin(Answer, thin = 200),x,f_0truth_at_x))
}
RatioObject=RatioObject*0.98
list_samples=list()
l=1
M=rep(0,8);
M[1]=0.06
for (y in 2:8){
M[y]=M[y-1]*RatioObject[y-1]
}
Count2=rep(0,8);
for (n in c(50,100,150,200,250,300,350,400)){
  
  Output=SimulationBM(5,1,1,n,0.2,0)
  Samples=Output[[1]]$trace;
  #f_theta=rep(0,nrows(Samples))
  for (h in 1:nrow(Samples)){
    distance=sqrt(mean((fourier_eval(as.vector(Samples[h,]),Output[[2]])-Output[[3]])^2))
    Count2[l]=Count2[l]+ifelse(distance<=M[l]*sqrt(5*log(n)/n),1,0)
  }
  print(Count2)
  l=l+1
}

n=100
l=1
M=0.08
s=1
g=1
j=1
ans_matrix=list()
for (t in 1:3){
  ans_matrix[[t]]=matrix(0,nrow=3,ncol=4)
}
coeffList<<-list(c(0.6,0.5),c(0.6,0.5,0.4,0.3),c(c(0.6,0.5,0.4,0.3,0.2)))
for (k in c(2,5,6)){
  g=1
  j=1
  for (alpha in c(0.5,1,1.5)){
    j=1
    for (tao in c(0.1,0.3,0.6,1)) {
      coeff <<- coeffList[[s]]
      Output=SimulationBM(k,alpha,tao,n,0.1, s)
      Samples=Output[[1]]$trace;
      print(Samples[1:2,])
      for (h in 1:nrow(Samples)){
        distance=sqrt(mean((fourier_eval(as.vector(Samples[h,]),Output[[2]])-Output[[3]])^2))
        ans_matrix[[s]][g,j]=ans_matrix[[s]][g,j]+ifelse(distance<=M*sqrt(k*log(n)/n),1,0)
      }
      print(ans_matrix)
      j=j+1
    }
    g=g+1
  }
  s=s+1
  
}




#
f_0truth=list()
f_0truth[[1]]<-function(x){
  coeff <- runif(5,0,3);
  res <- fourier_eval(coeff,x)
  return(res);
}
f_0truth[[2]]<-function(x){
  coeff <- runif(25,0,3);
  res <- fourier_eval(coeff,x)
  return(res);
}
f_0truth[[3]]<-function(x){
  coeff <- runif(50,0,3);
  res <- fourier_eval(coef,x)
  return(res);
}
MObject=rep(0,8)
RatioObject=rep(0,7)
s=1
for (jj in c(50,100,150,200,250,300,350,400)){
MObject[s]=(sqrt(5*log(jj)/jj))
s=s+1
}
for (d in 1:(length(MObject)-1)){
RatioObject[d]=MObject[d]/MObject[d+1];
}

for (n in c(50,100,150,200,250,300,350,400)){
print(M[l]*sqrt(5*log(n)/n))
l=l+1
}