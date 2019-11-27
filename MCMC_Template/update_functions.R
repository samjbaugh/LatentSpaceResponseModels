require('scatterpie')

sigmoid <- function(x)
{
  return(1/(1+exp(-x)))
}

euc_dist <- function(z,w)
{
  return(as.matrix(pdist(z,w)))
}

euc_dist_ordinal <- function(z,w)
{
  distout=pdist(z,w)
  return(array(as.matrix(distout),dim=c(nz,nw,ntau)))
}

source('plot_functions.R')

global_update<-function(name,axis,newvalue)
{
  myvariable=get(name)
  myvariable[[axis]]=newvalue
  assign(name,myvariable,envir=.GlobalEnv)
}

update_vector<-function(varname)
{
  current_vector=current_values[[varname]]
  n1=dim(current_vector)[1]
  n2=dim(current_vector)[2]
  proposal_sig=proposal_sigs[[varname]]
  proposed_vector=matrix(rnorm(n=n1*n2,mean=current_vector,sd=proposal_sig),n1,n2)
  
  likelihood_fun<-function(x) {likelihood_funs[[varname]](x)}
  prior_fun<-function(x) {prior_funs[[varname]](x)}
  
  logdiff_likelihood=likelihood_fun(proposed_vector)-likelihood_fun(current_vector)
  logdiff_prior=0 #prior_fun(proposed_vector)-prior_fun(current_vector)
  
  alpha=exp(logdiff_likelihood+logdiff_prior)
  randnum=runif(n1) #one for each k or i
  accepts=randnum<alpha
  
  retval=matrix(ifelse(rep(accepts,n2),proposed_vector,current_vector),n1,n2)
  if(any(is.na(retval)))
  {
    save(list=ls(),file='temp_save_delete.Rdat')
    error(e)
  }  
  return(list("newvalue"=retval,"alpha"=alpha,"accepts"=accepts)) 
}

update_is_spike<-function(varname="")
{
  likelihood_spike0=likelihood_logscale(current_values$logscale)
  likelihood_spike1=likelihood_logscale(-Inf)
  return(sample(c(0,1),1,p=c(likelihood_spike0,likelihood_spike1)))
}

update_K <- function(varname="")
{
  lambdas=current_values$lambda
  mus=current_values$mu
  sigmas=current_values$sigma
  
  probs_w=apply(current_values$w,1,function(xx) dnorm(xx[1],mus[,1],sigmas)*dnorm(xx[2],mus[,2],sigmas)*lambdas)
  probs_z=apply(current_values$z,1,function(xx) dnorm(xx[1],mus[,1],sigmas)*dnorm(xx[2],mus[,2],sigmas)*lambdas)
  
  sample_k<-function(p)
  {
    p[is.na(p)]=0
    p[is.null(p)]=0
    p[is.infinite(p)]=1
    if(all(p<1e-16))
    {
      p=rep(1,ncluster)
    }
    return(sample(size=1,x=1:ncluster,replace=T,prob=p))
  }
  K_z=apply(probs_z,2,sample_k)
  K_w=apply(probs_w,2,sample_k)
  
  return(list("K_z"=K_z,"K_w"=K_w))
}

update_mu_clusters <- function()
{
  sigs=current_values$sigma
  gms=stored_vars$gm
  gmeans=stored_vars$gmeans
  
  return(cbind(rnorm(ncluster,mean=gms*gmeans[,1]/(gms+sigs^2/current_values$omega^2),sd=sqrt(sigs^2/(gms+sigs^2/current_values$omega^2))),
               rnorm(ncluster,mean=gms*gmeans[,2]/(gms+sigs^2/current_values$omega^2),sd=sqrt(sigs^2/(gms+sigs^2/current_values$omega^2))))) ###take square root of sigmasqured term?
}


update_sigma <- function(varname)
{
  current_vector<-current_values[[varname]]
  n=length(current_vector)
  
  newsigma=sqrt(rinvgamma(1,shape=hyperparameters$invgam_shape_all+n/2,rate=hyperparameters$invgam_rate_all+sum(current_vector^2)/2))
  
  return(newsigma)
}


update_sigma_clusters <- function()
{
  gsd=stored_vars$gsd
  m=stored_vars$gm
  
  latent=1
  newsigmas=sqrt(rinvgamma(n=ncluster,shape=hyperparameters$invgam_shape_all+m*(latent+1)/2,rate=hyperparameters$invgam_rate_all+gsd^2/2))
  
  return(newsigmas)
}

update_lambda <- function()
{
 m=stored_vars$gm
 return(rdirichlet(1,m+hyperparameters$nu))
}

update_omega <- function()
{
  m=dim(current_values$mu)[1]
  new_omega=sqrt(rinvgamma(1,shape=hyperparameters$invgam_shape_all+m,rate=hyperparameters$invgam_rate_all+sum(current_values$mu^2)/4))
  return(new_omega)
}

update_stored_vars<-function()
{
  #store values for new Ks:
  gms=rep(NA,ncluster)
  gmeans=matrix(NA,ncluster,2)
  gsd=rep(NA,ncluster)
  for(ii in 1:ncluster)
  {
    temp=rbind(current_values$w[current_values$K_w==ii,],current_values$z[current_values$K_z==ii,])
    gms[ii]=dim(temp)[1]
    gmeans[ii,]=apply(temp,2,mean)
    gsd[ii]=sqrt(sum((temp-current_values$mu[ii,])^2)/2)
    if(gms[ii]==0)
    {
      gmeans[ii,]=c(0,0)
      gsd[ii]=0
    }
  }
  return(list('gm'=gms,'gmeans'=gmeans,'gsd'=gsd))
}


rdirichlet<-function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}
