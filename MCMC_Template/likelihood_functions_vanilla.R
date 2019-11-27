likelihood_funs=list()

likelihood_z<-function(given_z)
{
  myscale=c(exp(current_values$logscale))*(1-current_values$is_spike)
  wz_dist=myscale*euc_dist(given_z,current_values$w)
  bt_mat=outer(c(current_values$theta),c(current_values$beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(rowSums(retval))
}
likelihood_funs$z=likelihood_z

likelihood_w<-function(given_w)
{
  myscale=c(exp(current_values$logscale))*(1-current_values$is_spike)
  wz_dist=myscale*euc_dist(current_values[['z']],given_w)
  bt_mat=outer(c(current_values$theta),c(current_values$beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(colSums(retval))
}
likelihood_funs$w=likelihood_w

likelihood_theta<-function(given_theta)
{
  myscale=c(exp(current_values$logscale))*(1-current_values$is_spike)
  wz_dist=myscale*euc_dist(current_values$z,current_values$w)
  bt_mat=outer(c(given_theta),c(current_values$beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(rowSums(retval))
}
likelihood_funs$theta=likelihood_theta

likelihood_beta<-function(given_beta)
{
  myscale=c(exp(current_values$logscale))*(1-current_values$is_spike)
  wz_dist=myscale*euc_dist(current_values$z,current_values$w)
  bt_mat=outer(c(current_values$theta),c(given_beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(colSums(retval))
}
likelihood_funs$beta=likelihood_beta


likelihood_logscale<-function(given_logscale)
{
  if(current_values$is_spike)
  {
    return(-Inf)
  }
  myscale=c(exp(given_logscale))
  wz_dist=myscale*euc_dist(current_values$z,current_values$w)
  bt_mat=outer(c(current_values$theta),c(current_values$beta),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(sum(retval))
}
likelihood_funs$logscale=likelihood_logscale

calculate_full_likelihood<-function(stores,kk)
{
  wz_dist=c(exp(stores$logscale[[kk]]))*euc_dist(stores$z[[kk]],stores$w[[kk]])
  bt_mat=outer(c(stores$theta[[kk]]),c(stores$beta[[kk]]),'+')
  retval=sigmoid((2*X-1)*(bt_mat-wz_dist))
  return(sum(retval))
}


regular_prior<-function(given_vector,varname)
{
  return(-(given_vector^2)/(2*current_values[[paste('sigma',varname,sep='_')]]^2))
}

prior_funs=list()
prior_funs$beta=function(x) regular_prior(x,'beta')
prior_funs$theta=function(x) regular_prior(x,'theta')
prior_funs$w=function(x) regular_prior(x,'w')
prior_funs$z=function(x) regular_prior(x,'z')
prior_funs$logscale=function(x) regular_prior(x,'logscale')
