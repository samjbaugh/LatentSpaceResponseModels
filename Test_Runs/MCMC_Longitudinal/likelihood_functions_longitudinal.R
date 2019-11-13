likelihood_funs=list()

likelihood_z1<-function(given_z)
{
  wz_dist=euc_dist(given_z,current_values$w)
  bt_mat=outer(c(current_values$theta1),c(current_values$beta),'+')
  retval=sigmoid((2*X1-1)*(bt_mat-wz_dist))
  return(rowSums(log(retval)))
}
likelihood_funs$z1=likelihood_z1

likelihood_z2<-function(given_z)
{
  wz_dist=euc_dist(given_z,current_values$w)
  bt_mat=outer(c(current_values$theta2),c(current_values$beta),'+')
  retval=sigmoid((2*X2-1)*(bt_mat-wz_dist))
  return(rowSums(log(retval)))
}
likelihood_funs$z2=likelihood_z2

likelihood_w<-function(given_w)
{
  wz_dist=euc_dist(current_values$z1,given_w)
  bt_mat=outer(c(current_values$theta1),c(current_values$beta),'+')
  retval1=sigmoid((2*X1-1)*(bt_mat-wz_dist))
  
  wz_dist=euc_dist(current_values$z2,given_w)
  bt_mat=outer(c(current_values$theta2),c(current_values$beta),'+')
  retval2=sigmoid((2*X2-1)*(bt_mat-wz_dist))
  
  return(colSums(log(retval1+retval2)))
}
likelihood_funs$w=likelihood_w

likelihood_beta<-function(given_beta)
{
  wz_dist=euc_dist(current_values$z1,current_values$w)
  bt_mat=outer(c(current_values$theta1),c(given_beta),'+')
  retval1=sigmoid((2*X1-1)*(bt_mat-wz_dist))
  
  wz_dist=euc_dist(current_values$z2,current_values$w)
  bt_mat=outer(c(current_values$theta2),c(given_beta),'+')
  retval2=sigmoid((2*X2-1)*(bt_mat-wz_dist))
  
  return(colSums(log(retval1+retval2)))
}
likelihood_funs$beta=likelihood_beta

likelihood_theta1<-function(given_theta)
{
  wz_dist=euc_dist(current_values$z1,current_values$w)
  bt_mat=outer(c(given_theta),c(current_values$beta),'+')
  retval=sigmoid((2*X1-1)*(bt_mat-wz_dist))
  return(rowSums(log(retval)))
}
likelihood_funs$theta1=likelihood_theta1

likelihood_theta2<-function(given_theta)
{
  wz_dist=euc_dist(current_values$z2,current_values$w)
  bt_mat=outer(c(given_theta),c(current_values$beta),'+')
  retval=sigmoid((2*X2-1)*(bt_mat-wz_dist))
  return(rowSums(log(retval)))
}
likelihood_funs$theta2=likelihood_theta2

calculate_full_likelihood<-function(stores,kk)
{
  wz_dist=euc_dist(stores$z1[[kk]],stores$w[[kk]])
  bt_mat=outer(c(stores$theta1[[kk]]),c(stores$beta[[kk]]),'+')
  retval1=sigmoid((2*X1-1)*(bt_mat-wz_dist))
  
  wz_dist=euc_dist(stores$z2[[kk]],stores$w[[kk]])
  bt_mat=outer(c(stores$theta2[[kk]]),c(stores$beta[[kk]]),'+')
  retval2=sigmoid((2*X2-1)*(bt_mat-wz_dist))
  
  return(sum(retval1+retval2))
}

cluster_prior<-function(given_vector,varname,Kname)
{
  Ks=current_values[[Kname]]
  mus=current_values$mu[Ks,]
  sigmas=current_values$sigma[Ks]
  prior_val=-rowSums((given_vector-mus)^2)/(2*sigmas^2)
  return(prior_val)
}

regular_prior<-function(given_vector,varname)
{
  return(-(given_vector^2)/(2*current_values[[paste('sigma',varname,sep='_')]]^2))
}

latent_expand_covarmat<-function(covmat_base)
{
  covmat_expand=matrix(0,2*ndim,2*ndim)
  covmat_expand[1:ndim,1:ndim]=covmat_base
  covmat_expand[(ndim+1):(2*ndim),(ndim+1):(2*ndim)]=covmat_base
  x<-c(t(matrix(c(1:ndim,(ndim+1):(2*ndim)),ndim,ndim)))
  P=diag(length(x))[x,]
  covmat_expand=P%*%covmat_expand%*%t(P)
  return(covmat_expand)
}
latent_condense_covmat<-function(covmat_expand)
{
  x<-c(t(matrix(c(1:ndim,(ndim+1):(2*ndim)),ndim,ndim)))
  P=diag(length(x))[x,]
  permexp=P%*%covmat_expand%*%t(P)
  return(permexp[1:ndim,1:ndim])
}
wishart_prior<-function(covmat,latent=FALSE)
{
  condensefun=latent_condense_covmat
  outdense=dWishart(condensefun(covmat),df=wishart_df,Sigma=V)
  return(outdense)
}

prior_funs=list()
prior_funs$beta=function(x) regular_prior(x,'beta')
prior_funs$w=function(x) regular_prior(x,'w')
prior_funs$theta1=function(x) regular_prior(x,'theta1')
prior_funs$theta2=function(x) regular_prior(x,'theta2')
prior_funs$z1=function(x) cluster_prior(x,'z1','K_z')
prior_funs$z2=function(x) cluster_prior(x,'z2','K_z')

