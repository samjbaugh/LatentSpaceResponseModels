if(!file.exists('Run_configs')){dir.create('Run_configs')}
nu_config=100
rate_all_config=0.001
shape_all_config=0.001
omega_init_config=.1
sigma_theta_config=.1
sigma_mu_config=1
sigma_beta_config=1
sigma_w_config=.1
cluster_sigma_config=.1
config_number=1
sigma_tau_config=.1 #only needed for ordinal
config_filename=paste('Run_configs/config_',config_number,'.Rdat',sep='')
overwrite=T
if(!file.exists(config_filename) | overwrite)
{
  save(nu_config,rate_all_config,shape_all_config,sigma_w_config,
       omega_init_config,sigma_theta_config,sigma_mu_config,
       sigma_beta_config,cluster_sigma_config,sigma_tau_config,
       file=config_filename)
  print(paste(config_filename,'created'))
}else
{
  print('error: file exists')
}
