sigmoid <- function(x)
{
  return(1/(1+exp(-x)))
}

euc_dist <- function(z,w)
{
  return(as.matrix(pdist(z,w)))
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

plot_latent_cluster <- function(stored_parameters,kk,mytitle="",plot_pie=FALSE,burn_in=0,plot_ideology=F,save=F,save_filename="")
{
  z=data.frame(stored_parameters$z[[kk]]) #data.frame(matrix(stored_parameters$z[M,],nz,2))
  names(z)<-c('coord1','coord2')
  
  w=data.frame(stored_parameters$w[[kk]]) #data.frame(matrix(stored_parameters$w[M,],nw,2))
  names(w)<-c('coord1','coord2')
  
  z$K=gender
  
  mu=data.frame(stored_parameters$mu[[kk]])
  colnames(mu)<-c('coord1','coord2')
  w$name=sapply(1:nw,function(x) toString(x))
  mu$K=c('female','male')
  mu$r=stored_parameters$sigma[[kk]]
  
  p0<-ggplot()+geom_point(aes(x=coord1,y=coord2,col=K),data=z,cex=2)
  p0<-p0+xlab('coordinate 1')+ylab('coordinate 2')+ggtitle(paste('Latent space sample at k=',mytitle,sep=''))
  p0<-p0+geom_point(aes(x=coord1,y=coord2),pch=8,col='black',data=mu)
  p0<-p0+geom_circle(aes(x0=coord1,y0=coord2,r=r),data=mu)
  p0<-p0+geom_text(aes(x=coord1,y=coord2,label=name),hjust=2,vjust=2,data=w,fontface="bold")
  
  if(save)
  {
    png(save_filename)
    print(p0)
    dev.off()
  }else
  {
    print(p0)
  }
}

save_plot<-function(myp,name="myname"){
  png(filename=paste(name,".png",sep=""))
  print(myp)
  dev.off()
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
  
  logdiff_likelihood=likelihood_fun(proposed_vector)-likelihood_fun(current_vector) #done for numerical purposes
  logdiff_prior=prior_fun(proposed_vector)-prior_fun(current_vector)
  
  alpha=exp(logdiff_likelihood+logdiff_prior)
  randnum=runif(n1) #one for each k or i
  accepts=randnum<alpha
  
  retval=matrix(ifelse(rep(accepts,n2),proposed_vector,current_vector),n1,n2)
  
  return(list("newvalue"=retval,"alpha"=alpha,"accepts"=accepts)) 
}

# update_K <- function(varname="")
# {
#   lambdas=current_values$lambda
#   mus=current_values$mu
#   sigmas=current_values$sigma
#   
#   probs_w=apply(current_values$w,1,function(xx) dnorm(xx[1],mus[,1],sigmas)*dnorm(xx[2],mus[,2],sigmas)*lambdas)
#   probs_z=apply(current_values$z,1,function(xx) dnorm(xx[1],mus[,1],sigmas)*dnorm(xx[2],mus[,2],sigmas)*lambdas)
#   
#   sample_k<-function(p)
#   {
#     p[is.na(p)]=0
#     p[is.null(p)]=0
#     p[is.infinite(p)]=1
#     if(all(p<1e-16))
#     {
#       p=rep(1,ncluster)
#     }
#     return(sample(size=1,x=1:ncluster,replace=T,prob=p))
#   }
#   K_z=apply(probs_z,2,sample_k)
#   K_w=apply(probs_w,2,sample_k)
#   
#   return(list("K_z"=K_z,"K_w"=K_w))
# }

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
    # temp=rbind(current_values$w[current_values$K_w==ii,],current_values$z[current_values$K_z==ii,])
    temp=current_values$z[current_values$K_z==ii,]
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
