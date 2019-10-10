nu_config=2
rate_all_config=0.001
shape_all_config=0.001
sigma_tau_config=1
omega_config=sqrt(36)

stored_vars_init=list()
stored_vars_init[['n_z']]=nz
stored_vars_init[['n_w']]=nw
stored_vars_init[['ncluster']]=ncluster

sigma_theta_config=2
sigma_mu_config=2
sigma_beta_config=1
sigma_tau_config=1
cluster_sigma_config=1
config_number=1

if(!file.exists('Run_configs')) {dir.create(file.path('Run_configs'))}
config_filename=paste('Run_configs/config',config_number,'.Rdat',sep='')
overwrite=T
if(!file.exists(config_filename) | overwrite)
{
  save(nu_config,rate_all_config,shape_all_config,
       omega_config,sigma_theta_config,sigma_mu_config,
       sigma_beta_config,cluster_sigma_config,
       file=paste('Run_configs/config_',config_number,'.Rdat',sep=''))
}else
{
  print('error: file exists')
}
