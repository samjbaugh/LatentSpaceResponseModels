
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
  w_flattened=matrix(w,nw*ntau,2)
  distout=pdist(z,w_flattened)
  return(array(as.matrix(distout),dim=c(nz,nw,ntau)))
}

global_update<-function(name,axis,newvalue)
{
  myvariable=get(name)
  myvariable[[axis]]=newvalue
  assign(name,myvariable,envir=.GlobalEnv)
}


plot_latent <- function(stored_parameters,kk,title="")
{
  z=data.frame(stored_parameters$z[[kk]])
  w=data.frame(stored_parameters$w[[kk]])
  p0<-ggplot()
  p0<-p0+geom_point(aes(x=coord1,y=coord2,col="w"),data=z,cex=3)
  p0<-p0+geom_point(aes(x=coord1,y=coord2,col="z"),data=w,cex=3)
  p0<-p0+xlab('coordinate 1')+ylab('coordinate 2')+ggtitle(paste('Latent space sample at iteration',kk))
  print(p0)
}

plot_latent_ordinal_cluster <- function(stored_parameters,kk,save_tf=F,mytitle="",save_filename="")
{
  z=data.frame(stored_parameters$z[[kk]]) #data.frame(matrix(stored_parameters$z[M,],nz,2))
  w=data.frame(matrix(stored_parameters$w[[kk]],nw*ntau,2)) #data.frame(matrix(stored_parameters$w[M,],nw,2))
  names(z)<-c('coord1','coord2')
  names(w)<-c('coord1','coord2')
  
  temp_z=sapply(stored_parameters$K_z[[kk]],function(x) paste(toString(x)))
  z$K=as.factor(temp_z)
  
  temp_w=sapply(stored_parameters$K_w[[kk]],function(x) paste(toString(x)))
  w$K=as.factor(temp_w)
  w$wname=c(rep("w k=1",nw),rep("w k=2",nw),rep("w k=3",nw),rep("w k=4",nw))
  
  mu=data.frame(stored_parameters$mu[[kk]])
  colnames(mu)<-c('coord1','coord2')
  mu$wz="mu"
  mu$K=sapply(1:ncluster,function(x) paste(toString(x)))
  mu$r=stored_parameters$sigma[[kk]]
  
  p0<-ggplot()
  p0<-p0+geom_point(aes(x=coord1,y=coord2,col=K,pch='z'),data=z,cex=2)
  p0<-p0+geom_point(aes(x=coord1,y=coord2,col=K,pch=wname),data=w,cex=2)
  p0<-p0+geom_point(aes(x=coord1,y=coord2,col=K,pch='mu'),data=mu)
  p0<-p0+xlab('coordinate 1')+ylab('coordinate 2')+ggtitle(paste('Latent space sample at M=',mytitle,sep=''))
  p0<-p0+geom_circle(aes(x0=coord1,y0=coord2,col=K,r=r),data=mu)
  p0<-p0+scale_shape_manual(values=c("z"=16,"w k=1"=8,"w k=2"=9,"w k=3"=10,"w k=4"=11,"mu"=4))
  # p1=p0+geom_point(aes(x=current_values[['mu_z']][,1],y=current_values[['mu_z']][,2],col=c("1_z","2_z","3_z","4_z","5_z"),pch="cluster mu z"),cex=4)
  if(save_tf)
  {
    png(save_filename)
    print(p0)
    dev.off()
  }else
  {
    print(p0)
  }
}


update_vector<-function(varname)
{
  current_vector=current_values[[varname]]
  n1=dim(current_vector)[1]
  n2=dim(current_vector)[2]
  proposal_sig=proposal_sigs[[varname]]
  proposed_vector=matrix(rnorm(n=n1*n2,mean=current_vector,sd=proposal_sig),n1,n2)
  if(varname=='tau'){proposed_vector[,1]=0}
  
  likelihood_fun<-function(x) {likelihood_funs[[varname]](x)}
  prior_fun<-function(x) {prior_funs[[varname]](x)}
  
  logdiff_likelihood=likelihood_fun(proposed_vector)-likelihood_fun(current_vector) 
  logdiff_prior=prior_fun(proposed_vector)-prior_fun(current_vector)
  
  alpha=exp(logdiff_likelihood+logdiff_prior)
  randnum=runif(n1) #one for each k or i
  accepts=randnum<alpha
  
  retval=matrix(ifelse(rep(accepts,n2),proposed_vector,current_vector),n1,n2)
  
  return(list("newvalue"=retval,"alpha"=alpha,"accepts"=accepts)) 
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