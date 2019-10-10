initialize_sampler<-function(config_number,ordinal=F)
{
  if(!file.exists('Saved_output')) {dir.create('Saved_output')}
  if(!file.exists('Images')) {dir.create('Images')}
  
  load(paste('Run_configs/config_',config_number,'.Rdat',sep=''))
  
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
  
  K_z_init=sample(1:ncluster,nz,rep=T,prob=lambda_init)
  K_w_init=sample(1:ncluster,nw,rep=T,prob=lambda_init)
  
  mu_init=matrix(rnorm(2*ncluster,sd=sigma_mu_config),ncluster,2)
  sigma_init=rep(cluster_sigma_config,ncluster) #sqrt(rinvgamma(ncluster,shape=invgam_shape_latent,rate=invgam_rate_latent))   
  
  z_init=matrix(NA,nz,2)
  colnames(z_init)<-c('coord1','coord2')
  
  w_init=matrix(NA,nw,2)
  colnames(w_init)<-c('coord1','coord2')
  
  gms=rep(NA,ncluster)
  gmeans=matrix(NA,ncluster,2)
  gsd=rep(NA,ncluster)
  for(ii in 1:ncluster)
  {
    gmw=sum(K_w_init==ii)
    gmz=sum(K_z_init==ii)
    gm=gmw+gmz
    wii=cbind(rnorm(gmw,mean=mu_init[ii,1],sd=sigma_init[ii]),
              rnorm(gmw,mean=mu_init[ii,2],sd=sigma_init[ii]))
    zii=cbind(rnorm(gmz,mean=mu_init[ii,1],sd=sigma_init[ii]),
              rnorm(gmz,mean=mu_init[ii,2],sd=sigma_init[ii]))
    
    w_init[K_w_init==ii,]=wii
    z_init[K_z_init==ii,]=zii
    
    temp=rbind(wii,zii)
    gms[ii]=gm
    gmeans[ii,]=apply(temp,2,mean)
    gsd[ii]=sqrt(sum((temp-gmeans[ii,])^2)/2)
  }
  
  stored_vars_init=list()
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
                   "lambda"=lambda_init,
                   "K_z"=K_z_init,
                   "K_w"=K_w_init,
                   "omega"=omega_init)
  
  if(ordinal)
  {
    sigma_tau_init=sigma_tau_config
    tau_init=matrix(rnorm(nw*ntau,0,sigma_tau_init),nw,ntau)
    tau_init[,1]=0
    
    tau_mat=aperm(array(rep(tau_init,nz),dim=c(nw,ntau,nz)),perm=c(3,1,2))
    stored_vars_init$tau_mat=tau_mat
    
    init_values$tau=tau_init
    init_values$sigma_tau=sigma_tau_init
    
    assign("varname_list",c("z","w","theta","beta","tau"),envir=.GlobalEnv)
    assign("update_sigma_tf",list("z"=FALSE,"w"=FALSE,"theta"=TRUE,"beta"=FALSE,"tau"=TRUE),envir=.GlobalEnv)
  }else{
    assign("varname_list",c("z","w","theta","beta"),envir=.GlobalEnv)
    assign("update_sigma_tf",list("z"=FALSE,"w"=FALSE,"theta"=TRUE,"beta"=FALSE),envir=.GlobalEnv)
  }
  
  proposal_sigs=list()
  for(varname in varname_list)
  {
    proposal_sigs[[varname]]=2
  }
  
  return(list('init_values'=init_values,'stored_vars'=stored_vars_init,'hyperparameters'=hyperparameters,'proposal_sigs'=proposal_sigs))
}




