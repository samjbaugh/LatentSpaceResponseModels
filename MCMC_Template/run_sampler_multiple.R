
seedlist=c(6033,1833,6386,9986,9732,9643,1440,7670,6634,3071)

for(myseed in seedlist)
{
  print(myseed)
  run_mcmc_sampler(10000,myseed=myseed,config=1,plot_iter=100)
}
      