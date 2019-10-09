require('invgamma')
require('ggplot2')
require('pdist')
require('ggforce')
source('update_functions_cluster.R')
source('likelihood_functions.R')
source('init_clusters.R')
source('data_funs.R')

run_mcmc_sampler<-function(M,myseed,config_number,plot_iter=100,load_data=load_spelling_data)
{
  set.seed(myseed)
  #read data
  data_out=load_data()
  rawdata=data_out[['data']]
  dataname=data_out[['dataname']]
  gender<<-factor(rawdata$male)
  X<<-rawdata[,2:5]
  
  nz<<-dim(X)[1]
  nw<<-dim(X)[2]
  
  ndim<<-2
  ncluster<<-2
  
  #initalize
  init_out=initialize_cluster_sampler(config_number)
  current_values<<-init_out[['init_values']]
  stored_vars<<-init_out[['stored_vars']]
  hyperparameters<<-init_out[['hyperparameters']]
  
  varname_list=c("z","w","theta","beta") 
  update_sigma_tf=list("z"=FALSE,"w"=FALSE,"theta"=TRUE,"beta"=FALSE)
  
  plot_dirname=file.path(paste('Images/plots_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
  if(!file.exists(plot_dirname)) {dir.create(plot_dirname)}
  save_filename=file.path(paste('Saved_output/saved_output_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
  
  #initilize storage
  stored_parameters=list()
  
  for(varname in names(current_values))
  {
    stored_parameters[[varname]]=list()
    stored_parameters[[varname]][[1]]=current_values[[varname]]
  }
  
  stored_likelihoods=rep(NA,M)
  stored_likelihoods[1]=calculate_full_likelihood(stored_parameters,1)
  
  acceptance_rates=list()
  for(varname in c('z','w','beta','theta'))
  {
    acceptance_rates[[varname]]=rep(NA,M)
    acceptance_rates[[varname]][[1]]=1
  }
  global_vars=list()
  initial_proposal_densities=list('z'=2,'w'=2,'beta'=2,'theta'=2)
  for(varname in varname_list)
  {
    stored_parameters[[varname]][[1]]=current_values[[varname]]
    acceptance_rates[[varname]]=rep(NA,M)
    signame=paste('sigma',varname,sep='_')
    stored_parameters[[signame]][[1]]=current_values[[signame]]
    global_vars[[paste('proposal_sig_',varname,sep='')]]=initial_proposal_densities[[varname]]
  }
  global_vars<<-global_vars
  
  ###sample
  pb <- txtProgressBar(max = M, style = 3)
  plot_latent_cluster(stored_parameters,1,mytitle=toString(1),plot_pie=F,save=T,save_filename=paste(plot_dirname,'/iteration_',1,'.png',sep=''))
  
  # load('recent_save.Rdat')
  # start_index=length(stored_parameters$theta)
  plot_iter=100
  store_iter=1
  batch_size=500
  
  start_index=2
  runtimes=system.time({
    for(jj in start_index:M)
    {
      
      store=jj%%store_iter==0
      storej=floor(jj/store_iter)
      
      setTxtProgressBar(pb, jj)
      
      for(varname in varname_list)
      {
        out=update_vector(varname)
        
        current_values[[varname]]=out$newvalue
        if(store) {stored_parameters[[varname]][[storej]]=out$newvalue}
        
        acceptance_rates[[varname]][jj]=mean(out$accepts)
        
        if(update_sigma_tf[[varname]])
        {
          newsigma=update_sigma(varname)
          sig_varname=paste("sigma",varname,sep="_")
          if(store) {stored_parameters[[sig_varname]][[storej]]=newsigma}
          current_values[[sig_varname]]=newsigma
        }
      }
      newK=update_K()
      current_values$K_z=newK$K_z
      current_values$K_w=newK$K_w
      
      #store values for new Ks:
      stored_vars=update_stored_vars()
      
      newmu=update_mu_clusters()
      current_values$mu=newmu
      
      newsigma=update_sigma_clusters()
      current_values$sigma=newsigma
      
      newlambda=update_lambda()
      current_values$lambda=newlambda
      
      newomega=update_omega()
      current_values$omega=newomega
      
      stored_likelihoods[jj]=calculate_full_likelihood(stored_parameters,jj)
      
      if(store) {      
        stored_parameters$K_z[[storej]]=newK$K_z
        stored_parameters$K_w[[storej]]=newK$K_w
        stored_parameters$mu[[storej]]=newmu
        stored_parameters$sigma[[storej]]=newsigma
        stored_parameters$lambda[[storej]]=newlambda
        stored_parameters$omega[[storej]]=newomega
      }
      
      mi=max(1,jj-batch_size)
      for(varname in varname_list)
      {
        #rosenthall and roberts algorithm
        delta=ifelse(mean(acceptance_rates[[varname]][mi:jj],na.rm=T)>.44,min(1.01,1/sqrt(jj)),-min(1.01,1/sqrt(jj)))
        global_vars[[paste('proposal_sig_',varname,sep='')]]=global_vars[[paste('proposal_sig_',varname,sep='')]]*exp(delta)
      }
      
      if(jj%%plot_iter==0)
      {
        plot_latent_cluster(stored_parameters,storej,mytitle=toString(jj),plot_pie=F,save=T,save_filename=paste(plot_dirname,'/iteration_',jj,'.png',sep=''))
        save(stored_parameters,stored_likelihoods,global_vars,acceptance_rates,file=save_filename)   
      }
    }})
}

M=10000
myseed=1111
config_number=1
plot_iter=10
load_data=load_spelling_data
run_mcmc_sampler(M,myseed,config_number,plot_iter=plot_iter,load_data=load_data)
