sigmoid <- function(x)
{
  return(1/(1+exp(-x)))
}

euc_dist <- function(z,w)
{
  return(as.matrix(pdist(z,w)))
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
  z_indices=stored_parameters$K_z[[kk]] %in% c(1:3)
  z=z[z_indices,]
  
  w=data.frame(stored_parameters$w[[kk]]) #data.frame(matrix(stored_parameters$w[M,],nw,2))
  names(w)<-c('coord1','coord2')
  nz=dim(z)[1]
  nw=dim(w)[1]
  ncluster=2 #max(c(stored_parameters$K_z[[kk]],stored_parameters$K_w[[kk]]))
  
  temp_z=sapply(stored_parameters$K_z[[kk]][z_indices],function(x) paste(toString(x)))
  z$K=as.factor(temp_z)
  if(plot_ideology)
  {
    z$ideology=ideology
  }
  z$name=""

  temp_w=sapply(stored_parameters$K_w[[kk]],function(x) paste(toString(x)))
  w$K=as.factor(temp_w)
  w$name=sapply(1:nw,toString)
  
  total=length(stored_parameters$K_w)
  if(plot_pie)
  {
    clust_dist <- function(wnum) {return(c("1"=mean(sapply(stored_parameters$K_w[burn_in:total],function(x) x[wnum]==1)),"2"=mean(sapply(stored_parameters$K_w[burn_in:total],function(x) x[wnum]==2)),"3"=mean(sapply(stored_parameters$K_w[burn_in:total],function(x) x[wnum]==3))))}
    clust_totals <- sapply(1:nw,clust_dist)
    w$cluster1=clust_totals[1,]
    w$cluster2=clust_totals[2,]
    w$cluster3=clust_totals[3,]
  }

  mu=data.frame(stored_parameters$mu[[kk]])
  colnames(mu)<-c('coord1','coord2')
  mu$K=sapply(1:ncluster,function(x) paste(toString(x)))
  mu$r=stored_parameters$sigma[[kk]]

  if(plot_ideology)
  {
    p0<-ggplot()+geom_point(aes(x=coord1,y=coord2,col=K,pch=ideology),z,cex=2) #with ideology
    p0<-p0+scale_shape_manual(values=c("Liberal"=5,"Conservative"=16,"Moderate"=3))
  }else
  {
    p0<-ggplot()+geom_point(aes(x=coord1,y=coord2,pch=K,col=gender),z,cex=2)
  }
  p0<-p0+xlab('coordinate 1')+ylab('coordinate 2')+ggtitle(paste('Latent space sample at k=',mytitle,sep=''))
  if(ncluster<5)
  {
    p0<-p0+geom_point(aes(x=coord1,y=coord2),data=mu)
    p0<-p0+geom_circle(aes(x0=coord1,y0=coord2,r=r),data=mu)
  }
  if(plot_pie)
  {p0<-p0+geom_scatterpie(aes(x=coord1,y=coord2),data=w,cols=c("cluster1","cluster2","cluster3"))}
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
  ns=dim(current_vector) #stored_values[[paste('n',varname,sep='')]]
  n1=ns[1]
  n2=ns[2]
  proposal_sig=global_vars[[paste('proposal_sig',varname,sep='_')]]
  proposed_vector=matrix(rnorm(n=n1*n2,mean=current_vector,sd=proposal_sig),n1,n2)
  
  likelihood_fun<-function(x) {likelihood_funs[[varname]](x)}
  prior_fun<-function(x) {prior_funs[[varname]](x)}
  
  log_prob_proposed=likelihood_fun(proposed_vector)
  log_prob_current=likelihood_fun(current_vector)
  
  sumfun=switch(varname,'z'=rowSums,'theta'=rowSums,'w'=colSums,'beta'=colSums)
  
  logdiff_likelihood=sumfun(log_prob_proposed)-sumfun(log_prob_current) #done for numerical purposes
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

load_spelling_data<-function()
{
  return(list('data'=read.table('../Data/spelling.dat'),'dataname'='spelling'))
}

rdirichlet<-function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}
