require('invgamma')
require('ggplot2')
require('pdist')
require('ggforce')
source('update_functions_cluster.R')
source('likelihood_functions.R')
source('init_clusters.R')
source('data_funs.R')

run_mcmc_sampler<-function(M,myseed,config_number,plot_iter=1000,load_data=load_spelling_data)
{
  data_out=load_data()
  X<<-data_out[['data']]
  dataname<<-data_out[['dataname']]
  gender<<-data_out[['gender']]
  
  nz<<-dim(X)[1]
  nw<<-dim(X)[2]
  
  ndim<<-2
  ncluster<<-2
  
  store_iter=1
  batch_size=2000
  
  plot_dirname=file.path(paste('Images/plots_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
  if(!file.exists(plot_dirname)) {dir.create(plot_dirname)}
  save_filename=file.path(paste('Saved_output/saved_output_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
  
  varname_list<<-c("z","w","theta","beta") 
  
  overwrite=T
  if(file.exists(save_filename) & !overwrite)
  {
    load(save_filename,verb=T)
    start_index=length(stored_parameters$z)
    print(paste('Continuing from previous run at iteration:',start_index))
    
    set.seed(current_seed)
    assign("current_values",current_values,envir = .GlobalEnv)
    assign("stored_vars",stored_vars,envir = .GlobalEnv)
    assign("hyperparameters",hyperparameters,envir = .GlobalEnv)
    assign("proposal_sigs",proposal_sigs,envir = .GlobalEnv)
  }else
  {
    print('Starting from scratch.')
    #initalize
    set.seed(myseed)
    init_out=initialize_cluster_sampler(config_number)
    assign("current_values",init_out[['init_values']],envir = .GlobalEnv)
    assign("stored_vars",init_out[['stored_vars']],envir = .GlobalEnv)
    assign("hyperparameters",init_out[['hyperparameters']],envir = .GlobalEnv)
    assign("proposal_sigs",init_out[['proposal_sigs']],envir = .GlobalEnv)
    
    acceptance_rates=list()
    for(varname in varname_list)
    {
      acceptance_rates[[varname]]=rep(NA,M)
    }
    
    update_sigma_tf=list("z"=FALSE,"w"=TRUE,"theta"=TRUE,"beta"=FALSE)
    
    stored_parameters=list()
    for(varname in names(current_values))
    {
      stored_parameters[[varname]]=list()
      stored_parameters[[varname]][[1]]=current_values[[varname]]
    }
    
    stored_likelihoods=rep(NA,M)
    stored_likelihoods[1]=calculate_full_likelihood(stored_parameters,1)
    
    ###sample
    plot_latent_cluster(stored_parameters,1,mytitle='Initial Configuration',plot_pie=F,save=T,save_filename=paste(plot_dirname,'/initial_configuration.png',sep=''))
    start_index=2
  }
  
  pb <- txtProgressBar(max = M, style = 3)
  
  for(jj in start_index:M)
  {
    current_seed=runif(1)*1e9 #for continuation purposes
    set.seed(current_seed)
    store=jj%%store_iter==0
    storej=floor(jj/store_iter)
    
    setTxtProgressBar(pb, jj)
    
    for(varname in varname_list)
    {
      out=update_vector(varname)
      
      global_update("current_values",varname,out$newvalue)
      
      if(store) {stored_parameters[[varname]][[storej]]=out$newvalue}
      
      acceptance_rates[[varname]][jj]=mean(out$accepts)
      
      if(update_sigma_tf[[varname]])
      {
        newsigma=update_sigma(varname)
        sig_varname=paste("sigma",varname,sep="_")
        if(store) {stored_parameters[[sig_varname]][[storej]]=newsigma}
        global_update("current_values",sig_varname,newsigma)
      }
    }
    assign("stored_vars",update_stored_vars(),envir = .GlobalEnv)
    
    # newK=update_K()
    # global_update("current_values","K_w",new$K_w)
    # global_update("current_values","K_z",new$K_z)
    
    #store values for new Ks:
    global_update("current_values","mu",update_mu_clusters())
    global_update("current_values","sigma",update_sigma_clusters())
    global_update("current_values","lambda",update_lambda())
    global_update("current_values","omega",update_omega())
    
    stored_likelihoods[jj]=calculate_full_likelihood(stored_parameters,jj)
    
    if(store) {      
      stored_parameters$K_z[[storej]]=current_values$K_z
      stored_parameters$mu[[storej]]=current_values$mu
      stored_parameters$sigma[[storej]]=current_values$sigma
      stored_parameters$lambda[[storej]]=current_values$lambda
      stored_parameters$omega[[storej]]=current_values$omega
    }
    
    
    mi=max(1,jj-batch_size)
    for(varname in varname_list)
    {
      #rosenthall and roberts algorithm
      delta=ifelse(mean(acceptance_rates[[varname]][1:jj],na.rm=T)>.44,min(.01,1/sqrt(jj)),-min(.01,1/sqrt(jj)))
      proposal_sigs[[varname]]=proposal_sigs[[varname]]*exp(delta)
    }
    assign("proposal_sigs",proposal_sigs,envir = .GlobalEnv)
    
    if(jj%%plot_iter==0)
    {
      # print('acceptance rates:')
      # print(sapply(acceptance_rates,function(x) mean(x,na.rm=T)))
      # print(data.frame(proposal_sigs))
      plot_latent_cluster(stored_parameters,storej,mytitle=toString(jj),save=T,save_filename=paste(plot_dirname,'/iteration_',jj,'.png',sep=''))
      save(stored_parameters,current_values,stored_likelihoods,stored_vars,hyperparameters,proposal_sigs,acceptance_rates,current_seed,file=save_filename)   
    }
  }
}

# M=100000
# myseed=2222
# config_number=1
# plot_iter=1000
# load_data=load_spelling_data
# run_mcmc_sampler(M,myseed,config_number,plot_iter=plot_iter,load_data=load_data)
# 
# save_filename=file.path(paste('Saved_output/saved_output_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
# load(save_filename,verb=T)
# print(sapply(acceptance_rates,function(x) mean(x,na.rm=T)))
