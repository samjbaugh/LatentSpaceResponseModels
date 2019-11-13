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
  if(ordinal)
  {
    K_w_init=sample(1:ncluster,nw*ntau,rep=T,prob=lambda_init)
  }else{
    K_w_init=sample(1:ncluster,nw,rep=T,prob=lambda_init)
  }
  
  mu_init=matrix(rnorm(2*ncluster,sd=sigma_mu_config),ncluster,2)
  sigma_init=rep(cluster_sigma_config,ncluster) #sqrt(rinvgamma(ncluster,shape=invgam_shape_latent,rate=invgam_rate_latent))   
  
  z_init=matrix(NA,nz,2)
  colnames(z_init)<-c('coord1','coord2')
  
  if(ordinal)
  {
    w_init=matrix(NA,nw*ntau,2)
  }
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
  if(ordinal)
  {
    sigma_beta=0
    beta_init=matrix(0,nw,1)
  }else{
    sigma_beta=sigma_beta_config
    beta_init=matrix(rnorm(nw,0,sigma_beta),nw,1)
  }

  
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
    
    assign("varname_list",c("z","w","theta","tau"),envir=.GlobalEnv)
    assign("update_sigma_tf",list("z"=FALSE,"w"=FALSE,"theta"=TRUE,"tau"=TRUE),envir=.GlobalEnv)
  }
  vanilla=!ordinal
  if(vanilla)
  {
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

########################################
##########LONGITUDINAL VERSION##########
########################################

initialize_sampler_longitudinal<-function(config_number,ordinal=F)
{
  if(!file.exists('Saved_output')) {dir.create('Saved_output')}
  if(!file.exists('Images')) {dir.create('Images')}
  
  load(paste('Run_configs/config_',config_number,'.Rdat',sep=''))
  
  hyperparameters=list()
  init_values=list()
  stored_vars_init=list()
  
  ######init hyperparameters######
  hyperparameters$nu=nu_config
  hyperparameters$invgam_shape_all=shape_all_config
  hyperparameters$invgam_rate_all=rate_all_config
  
  ######init cluser params######
  omega_init=omega_init_config
  init_values[['omega']]=omega_init
  lambda_init=rdirichlet(1,rep(hyperparameters$nu,ncluster))
  init_values[['lambda']]=lambda_init
  
  ######init theta######
  sigma_theta1_init=sigma_theta_config #sqrt(rinvgamma(1,shape=invgam_shape_additive,rate=invgam_rate_additive))
  sigma_theta2_init=sigma_theta_config #sqrt(rinvgamma(1,shape=invgam_shape_additive,rate=invgam_rate_additive))
  theta1_init=matrix(rnorm(nz,0,sigma_theta1_init),nz,1)
  theta2_init=matrix(rnorm(nz,0,sigma_theta2_init),nz,1)
  init_values[['sigma_theta1']]=sigma_theta1_init
  init_values[['theta1']]=theta1_init
  init_values[['sigma_theta2']]=sigma_theta2_init
  init_values[['theta2']]=theta2_init
  
  ######init z######
  K_z_init=sample(1:ncluster,nz,rep=T,prob=lambda_init)
  mu_init=matrix(rnorm(2*ncluster,sd=sigma_mu_config),ncluster,2)
  sigma_init=rep(cluster_sigma_config,ncluster) #sqrt(rinvgamma(ncluster,shape=invgam_shape_latent,rate=invgam_rate_latent))   
  init_values[['K_z']]=K_z_init
  init_values[['mu']]=mu_init
  init_values[['sigma']]=sigma_init
  
  z1_init=matrix(NA,nz,2)
  colnames(z1_init)<-c('coord1','coord2')
  z2_init=matrix(NA,nz,2)
  colnames(z2_init)<-c('coord1','coord2')

  gms=rep(NA,ncluster)
  gmeans=matrix(NA,ncluster,2)
  gsd=rep(NA,ncluster)
  for(ii in 1:ncluster)
  {
    gmz=sum(K_z_init==ii)
    gm=gmz
    #maybe make multivariate sample?
    z1_ii=cbind(rnorm(gmz,mean=mu_init[ii,1],sd=sigma_init[ii]),
              rnorm(gmz,mean=mu_init[ii,2],sd=sigma_init[ii]))
    z2_ii=cbind(rnorm(gmz,mean=mu_init[ii,1],sd=sigma_init[ii]),
               rnorm(gmz,mean=mu_init[ii,2],sd=sigma_init[ii]))   
    z1_init[K_z_init==ii,]=z1_ii
    z2_init[K_z_init==ii,]=z2_ii
    
    temp=rbind(z1_ii,z2_ii)
    gms[ii]=gm
    gmeans[ii,]=apply(temp,2,mean)
    #change to empirical covariance matrix
    gsd[ii]=sqrt(sum((temp-gmeans[ii,])^2)/2)
  }
  init_values[['z1']]=z1_init
  init_values[['z2']]=z2_init
  
  stored_vars_init$gm=gms
  stored_vars_init$gmeans=gmeans
  stored_vars_init$gsd=gsd
  
  ######init w######
  sigma_w_init=sigma_w_config
  w_init=cbind(rnorm(nw,mean=0,sd=sigma_w_init),
               rnorm(nw,mean=0,sd=sigma_w_init)) 
  colnames(w_init)<-c('coord1','coord2')
  init_values[['w']]=w_init
  init_values[['sigma_w']]=w_init
  
  ######init beta######
  sigma_beta=sigma_beta_config
  beta_init=matrix(rnorm(nw,0,sigma_beta),nw,1)
  init_values[['beta']]=beta_init
  init_values[['sigma_beta']]=sigma_beta
  
  ######init varname_list######
  assign("varname_list",c("z1","z2","w","theta1","theta2","beta"),envir=.GlobalEnv)
  assign("update_sigma_tf",list("z1"=FALSE,"z2"=FALSE,"theta1"=TRUE,"theta2"=TRUE,"w"=F,"beta"=FALSE),envir=.GlobalEnv)
  
  proposal_sigs=list()
  for(varname in varname_list)
  {
    proposal_sigs[[varname]]=2
  }
  
  return(list('init_values'=init_values,'stored_vars'=stored_vars_init,'hyperparameters'=hyperparameters,'proposal_sigs'=proposal_sigs))
}




