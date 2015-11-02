library(fda)
library(mvtnorm)
library(MHadaptive)
sin_basis<-function(x,n){
  return(sqrt(2)*sin(2*pi*n*x));  
}
#Sin Basis Function

cos_basis<-function(x,n){
  return(sqrt(2)*cos(2*pi*n*x));   
}
#Cos Basis Function

f_0truth<-function(x){
  f_0=2+1.2*sin_basis(x,1)+1*cos_basis(x,1)+0.8*sin_basis(x,2)+0.7*cos_basis(x,2);
  return(f_0)
}
#Truth function we aim to be near to.

SimulationBM<-function(k,alpha,tao,n,sigma){
  #Simulation of the quantity of interest with inputs:
  # k the number of basis (would like to be greater than number of basis of the f_0)
  # alpha, tao and sigma control
  # n the number of data
  
  x<<-runif(n,0,1);
  # Generate the covariates (n of these)
  
  #Ignore
  #thetaPri=rep(0,k);
  #for (j in 1:k) {
  #theta[j]=rnorm(1,0,sqrt((tao^2)*j^(-2*alpha-1)));
  #}
  
  epsilon=rnorm(n,0,sigma^2);
  #Generate the errors with sigma 
  
  f_0truth_at_x <<-f_0truth(x);
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
  
  
  #N=1e3
  #p=200
  #d=5
  
  #priorf = function(mu) {dnorm(mu[1],1.5,0.55) * dnorm(mu[2],1.5,0.55) * dnorm(mu[3],1.5,0.55) * dnorm(mu[4],1.5,0.55)}
  #priorf = function(mu) {dmvnorm(mu, mean=c(1.5,1.5,1.5,1.5), sigma=0.55*diag(4))}
  
  #mu is 4-d vector
  
  
  #assuming data y has been generated
  #likelic = function(mu) {prod(rowSums(0.25*(sapply(mu,function(x){dnorm(y,x,0.55)}))))}
  #likelic = function(mu) {}
  
  #phi = function(t){t^2} #Later change this to 6 or 7
  #tseq = seq(from = 0, to = 1, len = p)
  
  
  #sequenceBayes = function(i) {
   # z = i
  #  fn = function(theta) {
   #   return(((Likelihood(y_data,x,epsilon,sigma,n,f_0truth_at_x,theta))^z)*PriorMult(theta,tao,alpha,k))
  #  }
  #  return(fn)
  #}
  
#  seqBayes = list()
  #sample from prior:
 # rinit = rep(0,5);
  #rinit1 = t(replicate(n=N,c(rnorm(1,-3,0.55),rnorm(1,0,0.55),rnorm(1,3,0.55),rnorm(1,6,0.55))))
#  for(i in 1:p) {
 #   seqBayes[[i]] = sequenceBayes(phi(tseq[i]))
#  }
  #return(seqBayes)
  #return(sequentialMC(N,d,rinit,seqBayes))
         
         #sequentialMC = function(N, d, rinit, gamma.list, T = N/2, nMH = 10) {
         
         
#return(MetropolisHastingC(1000,function(theta) {PosteriorKernel(theta,y_data,x,epsilon,sigma,n,f_0truth_at_x,tao,alpha,k)},0.001,c(2,1,1,1,1),5))
         #Generates 1000 samples from posterior distribution, with independent multivariate gaussian as proposal with 0.01 as variance across diagnoal.
         #Intialisation is at 0 k dimensioanl with thinning at 5.

Answer=(Metro_Hastings(function(theta) {PosteriorKernel(theta,y_data,x,epsilon,sigma,n,f_0truth_at_x,tao,alpha,k)},c(1.5,2,1.5,1.5,1.5),iterations = 500*10+500,burn_in = 500))
return(mcmc_thin(Answer, thin = 10))
}
list_samples=list()
l=1
M=c(2.6,2.65,2.7,2.75,2.8);
Count=rep(0,5)
for (n in c(100,200,300,400,500)){
  Samples=SimulationBM(5,1.5,1,n,0.1)$trace;
  #f_theta=rep(0,nrows(Samples))
  for (h in 1:nrow(Samples)){
  Count[l]=Count[l]+ifelse((mean(fourier_eval(as.vector(Samples[h,]),x)-f_0truth_at_x)<=M[l]*sqrt((5*log(n)/n))),1,0)
  }
  Count=Count/nrow(Samples)
  l=l+1
}




















MetropolisHastingC = function(TamplesN, densityInt, sigma, initial,thin) {
  SamplesN=thin*TamplesN;
  #Works out how much thinning
  OldPos=initial;
  OutputSamples= matrix(0,SamplesN,length(initial));
  # pb <- txtProgressBar(min = 0, max = SamplesN, style = 3);
  #Allocation
  Dimensions<-length(initial);
  for (i in 1:SamplesN)
  {
    Proposal<<- rmvnorm(1,as.vector(OldPos),diag(sigma,Dimensions));
    #as.vector(OldPos)+as.vector(runif(Dimensions,-sigma,sigma))
    #rmvnorm(1,as.vector(OldPos),sigma);
    ProposedValue<<-log(densityInt(as.vector(Proposal)));
    #Proposed Value Density
    OldValue<<-log(densityInt(as.vector(OldPos)));
    #Old Value Density
    if ((runif(1))<= as.brob(exp(ProposedValue-OldValue)))
      #MH rejection
    {
      OldPos=Proposal;
    }
    OutputSamples[i,]=OldPos
    #setTxtProgressBar(pb, i)
  }
  if (TamplesN==1)
  {
    FinalOutputSamples=OutputSamples[thin,]
  }
  else {
    #Size=dim(unique(OutputSamples))
    FinalOutputSamples=OutputSamples[seq(1,dim(OutputSamples)[1],thin),] }
  #Thinning
  #plot(FinalOutputSamples[,1],FinalOutputSamples[,2],type='l',col=2,xlab="X",ylab="Y",xlim=c(-2,2),ylim=c(-2,2))
  #points(FinalOutputSamples[,1],FinalOutputSamples[,2],pch=4)
  #cat("rejectionRate=",1-(Size[1]/SamplesN))
  #after=Sys.time()-now;
  #print(after)
  return(FinalOutputSamples)
  #close(pb)
}
