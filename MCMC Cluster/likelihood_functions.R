likelihood_funs=list()

likelihood_z<-function(given_z)
{
  wz_dist=euc_dist(given_z,current_values$w)
  bt_mat=outer(c(current_values$theta),c(current_values$beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(rowSums(retval))
}
likelihood_funs$z=likelihood_z

likelihood_w<-function(given_w)
{
  wz_dist=euc_dist(current_values[['z']],given_w)
  bt_mat=outer(c(current_values$theta),c(current_values$beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(colSums(retval))
}
likelihood_funs$w=likelihood_w

likelihood_theta<-function(given_theta)
{
  wz_dist=euc_dist(current_values$z,current_values$w)
  bt_mat=outer(c(given_theta),c(current_values$beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(rowSums(retval))
}
likelihood_funs$theta=likelihood_theta

likelihood_beta<-function(given_beta)
{
  wz_dist=euc_dist(current_values$z,current_values$w)
  bt_mat=outer(c(current_values$theta),c(given_beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(colSums(retval))
}
likelihood_funs$beta=likelihood_beta

calculate_full_likelihood<-function(stores,kk)
{
  wz_dist=euc_dist(stores$z[[kk]],stores$w[[kk]])
  bt_mat=outer(c(stores$theta[[kk]]),c(stores$beta[[kk]]),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
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
