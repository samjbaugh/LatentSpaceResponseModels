
seedlist=c(6033,1833,6386,9986,9732,9643,1440,7670,6634,3071,623,24234,12451,12467,12313)

for(myseed in seedlist)
{
  print(myseed)
  run_mcmc_sampler(10000,myseed=myseed,config=1,plot_iter=100)
}
      

seedlist=c(6033,1833,6386,9986,9732,9643,1440,7670,6634,3071,623,24234,12451,12467,12313)

for(myseed in seedlist)
{
  save_filename=file.path(paste('Saved_output/saved_output_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
  load(save_filename,verb=T)
}