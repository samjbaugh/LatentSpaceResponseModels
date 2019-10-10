nu_config=10
rate_all_config=0.001
shape_all_config=0.001
omega_init_config=3
sigma_theta_config=2
sigma_mu_config=2
sigma_beta_config=1
cluster_sigma_config=1
config_number=1
config_filename=paste('Run_configs/config',config_number,'.Rdat',sep='')
overwrite=F
if(!file.exists(config_filename) | overwrite)
{
  save(nu_config,rate_all_config,shape_all_config,
       omega_init_config,sigma_theta_config,sigma_mu_config,
       sigma_beta_config,cluster_sigma_config,
       file=config_filename)
}else
{
  print('error: file exists')
}
