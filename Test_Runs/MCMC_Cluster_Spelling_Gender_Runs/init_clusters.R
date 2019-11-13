initialize_cluster_sampler<-function(config_number)
{
  load(paste('Run_configs/config_',config_number,'.Rdat',sep=""))
  
  stored_vars_init=list()
  stored_vars_init$n_z=nz
  stored_vars_init$n_w=nw
  stored_vars_init$ncluster=ncluster
  
  #hyperparameters
  hyperparameters=list()
  hyperparameters$nu=nu_config
  hyperparameters$invgam_shape_all=shape_all_config
  hyperparameters$invgam_rate_all=rate_all_config
  
  #configurable initialization parameters
  omega_init=omega_init_config
  
  sigma_theta_init=sigma_theta_config #sqrt(rinvgamma(1,shape=invgam_shape_additive,rate=invgam_rate_additive))
  theta_init=matrix(rnorm(nz,0,sigma_theta_init),nz,1)
  
  #initialize z:
  lambda_init=rdirichlet(1,rep(hyperparameters$nu,ncluster))
  
  K_z_init=as.numeric(gender)
  mu_init=matrix(rnorm(2*ncluster,sd=sigma_mu_config),ncluster,2)
  sigma_init=rep(cluster_sigma_config,ncluster) #sqrt(rinvgamma(ncluster,shape=invgam_shape_latent,rate=invgam_rate_latent))   
  
  z_init=matrix(NA,nz,2)
  colnames(z_init)<-c('coord1','coord2')
  
  sigma_w_init=sigma_w_config
  w_init=cbind(rnorm(nw,mean=0,sd=sigma_w_init),
                   rnorm(nw,mean=0,sd=sigma_w_init)) 
  colnames(w_init)<-c('coord1','coord2')
  
  gms=rep(NA,ncluster)
  gmeans=matrix(NA,ncluster,2)
  gsd=rep(NA,ncluster)
  for(ii in 1:ncluster)
  {
    gmz=sum(K_z_init==ii)
    gm=gmz
    zii=cbind(rnorm(gmz,mean=mu_init[ii,1],sd=sigma_init[ii]),
              rnorm(gmz,mean=mu_init[ii,2],sd=sigma_init[ii]))
    
    z_init[K_z_init==ii,]=zii
    
    temp=zii
    gms[ii]=gm
    gmeans[ii,]=apply(temp,2,mean)
    gsd[ii]=sqrt(sum((temp-gmeans[ii,])^2)/2)
  }
  stored_vars_init$gm=gms
  stored_vars_init$gmeans=gmeans
  stored_vars_init$gsd=gsd
  
  #init beta
  sigma_beta=sigma_beta_config
  beta_init=matrix(rnorm(nw,0,sigma_beta),nw,1)
  
  init_values=list("z"=z_init,
                   "w"=w_init,
                   "theta"=theta_init,
                   "sigma_theta"=sigma_theta_init,
                   "beta"=beta_init,
                   "sigma_beta"=sigma_beta,
                   "mu"=mu_init,
                   "sigma"=sigma_init,
                   "sigma_w"=sigma_w_init,
                   "lambda"=lambda_init,
                   "K_z"=K_z_init,
                   # "K_w"=K_w_init,
                   "omega"=omega_init)
  
  proposal_sigs=list()
  for(varname in varname_list)
  {
    proposal_sigs[[varname]]=2
  }
  
  return(list('init_values'=init_values,'stored_vars'=stored_vars_init,'hyperparameters'=hyperparameters,'proposal_sigs'=proposal_sigs))
}




