likelihood_funs=list()

likelihood_z<-function(given_z)
{
  tau_mat=aperm(array(rep(current_values$tau,nz),dim=c(nw,ntau,nz)),perm=c(3,1,2))
  
  wz_dist=euc_dist_ordinal(given_z,current_values$w)
  bt_mat=array(rep(outer(c(current_values$theta),c(current_values$beta),'+'),ntau),dim=c(nz,nw,ntau))
  
  term_mat_unnorm=sigmoid(tau_mat+bt_mat-wz_dist)
  term_mat_norm=aperm(apply(term_mat_unnorm,c(1,2),function(r) r/sum(r)),c(2,3,1)) #normalize
  
  retval=sapply(1:nw,function(jj) 
    sapply(1:nz,function(ii) term_mat_norm[ii,jj,][X[ii,jj]-minord+1]))
  retval[retval==0]=1e-100
  return(rowSums(log(retval)))
}
likelihood_funs$z=likelihood_z

likelihood_w<-function(given_w)
{
  tau_mat=aperm(array(rep(current_values$tau,nz),dim=c(nw,ntau,nz)),perm=c(3,1,2))
  
  wz_dist=euc_dist_ordinal(current_values$z,given_w)
  bt_mat=array(rep(outer(c(current_values$theta),c(current_values$beta),'+'),ntau),dim=c(nz,nw,ntau))
  
  term_mat_unnorm=sigmoid(tau_mat+bt_mat-wz_dist)
  term_mat_norm=aperm(apply(term_mat_unnorm,c(1,2),function(r) r/sum(r)),c(2,3,1)) #normalize
  
  retval=sapply(1:nw,function(jj) 
    sapply(1:nz,function(ii) term_mat_norm[ii,jj,][X[ii,jj]-minord+1]))
  retval[retval==0]=1e-100
  return(colSums(log(retval)))
}
likelihood_funs$w=likelihood_w

likelihood_theta<-function(given_theta)
{
  tau_mat=aperm(array(rep(current_values$tau,nz),dim=c(nw,ntau,nz)),perm=c(3,1,2))
  
  wz_dist=euc_dist_ordinal(current_values$z,current_values$w)
  bt_mat=array(rep(outer(c(given_theta),c(current_values$beta),'+'),ntau),dim=c(nz,nw,ntau))
  
  term_mat_unnorm=sigmoid(tau_mat+bt_mat-wz_dist)
  term_mat_norm=aperm(apply(term_mat_unnorm,c(1,2),function(r) r/sum(r)),c(2,3,1)) #normalize
  
  retval=sapply(1:nw,function(jj) 
    sapply(1:nz,function(ii) term_mat_norm[ii,jj,][X[ii,jj]-minord+1]))
  retval[retval==0]=1e-100
  return(rowSums(log(retval)))
}
likelihood_funs$theta=likelihood_theta

likelihood_beta<-function(given_beta)
{
  tau_mat=aperm(array(rep(current_values$tau,nz),dim=c(nw,ntau,nz)),perm=c(3,1,2))
  
  wz_dist=euc_dist_ordinal(current_values$z,current_values$w)
  bt_mat=array(rep(outer(c(current_values$theta),c(given_beta),'+'),ntau),dim=c(nz,nw,ntau))
  
  term_mat_unnorm=sigmoid(tau_mat+bt_mat-wz_dist)
  term_mat_norm=aperm(apply(term_mat_unnorm,c(1,2),function(r) r/sum(r)),c(2,3,1)) #normalize
  
  retval=sapply(1:nw,function(jj) 
    sapply(1:nz,function(ii) term_mat_norm[ii,jj,][X[ii,jj]-minord+1]))
  retval[retval==0]=1e-100
  return(colSums(log(retval)))
}
likelihood_funs$beta=likelihood_beta

likelihood_tau<-function(given_tau)
{
  tau_mat=aperm(array(rep(given_tau,nz),dim=c(nw,ntau,nz)),perm=c(3,1,2))
  
  wz_dist=euc_dist_ordinal(current_values$z,current_values$w)
  bt_mat=array(rep(outer(c(current_values$theta),c(current_values$beta),'+'),ntau),dim=c(nz,nw,ntau))
  
  term_mat_unnorm=sigmoid(tau_mat+bt_mat-wz_dist)
  term_mat_norm=aperm(apply(term_mat_unnorm,c(1,2),function(r) r/sum(r)),c(2,3,1)) #normalize
  
  retval=sapply(1:ntau,function(k) {r=term_mat_norm[,,k];
                r*(X==k)+(1-(X==k))*(1-r)},simplify="array")
  retval[retval==0]=1e-100
  return(apply(log(retval),c(2,3),sum))
  ##end likelihood section
}
likelihood_funs$tau=likelihood_tau


calculate_full_likelihood<-function(stored_parameters,kk)
{
  tau_mat=aperm(array(rep(stored_parameters$tau[[kk]],nz),dim=c(nw,ntau,nz)),perm=c(3,1,2))
  
  wz_dist=euc_dist_ordinal(stored_parameters$z[[kk]],stored_parameters$w[[kk]])
  bt_mat=array(rep(outer(c(stored_parameters$theta[[kk]]),c(stored_parameters$beta[[kk]]),'+'),ntau),dim=c(nz,nw,ntau))
  
  term_mat_unnorm=sigmoid(tau_mat+bt_mat-wz_dist)
  term_mat_norm=aperm(apply(term_mat_unnorm,c(1,2),function(r) r/sum(r)),c(2,3,1)) #normalize
  
  retval=sapply(1:nw,function(jj) 
    sapply(1:nz,function(ii) term_mat_norm[ii,jj,][X[ii,jj]-minord+1]))
  return(sum(retval))
}

cluster_prior<-function(given_vector,varname)
{
  Ks=current_values[[paste('K',varname,sep='_')]]
  mus=current_values$mu[Ks,]
  sigmas=current_values$sigma[Ks]
  prior_val=-rowSums((given_vector-mus)^2)/(2*sigmas^2)
  return(prior_val)
}

regular_prior<-function(given_vector,varname)
{
  return(-(given_vector^2)/(2*current_values[[paste('sigma',varname,sep='_')]]^2))
}

prior_funs=list()
prior_funs$beta=function(x) regular_prior(x,'beta')
prior_funs$theta=function(x) regular_prior(x,'theta')
prior_funs$w=function(x) cluster_prior(x,'w')
prior_funs$z=function(x) cluster_prior(x,'z')
prior_funs$tau=function(x) regular_prior(x,'tau')
