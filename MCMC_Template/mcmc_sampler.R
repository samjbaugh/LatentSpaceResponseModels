require('invgamma')
require('ggplot2')
require('pdist')
require('ggforce')
require('progress')
source('update_functions.R')
source('init_sampler.R')
source('data_funs.R')

run_mcmc_sampler<-function(M,myseed,config_number,plot_iter=1000,
                           load_data,ordinal,save_fig=T,ncluster=2,
                           overwrite=F,return_output=T)
{
  if(ordinal){
    source('likelihood_functions_ordinal.R')
    plot_fun=plot_latent_ordinal_cluster
  }else if(ncluster>1){
    source('likelihood_functions_cluster.R')
    plot_fun=plot_latent_cluster
  }else{
    source('likelihood_functions_vanilla.R')
    plot_fun=plot_latent
  }
  
  load_data()
  
  ndim<<-2
  assign("ncluster",ncluster,envir = .GlobalEnv)
  
  store_iter=1
  batch_size=2000
  
  plot_dirname=file.path(paste('Images/plots_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
  if(!file.exists(plot_dirname)) {dir.create(plot_dirname)}
  save_filename=file.path(paste('Saved_output/saved_output_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
  
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
    init_out=initialize_sampler(config_number,ordinal=ordinal)
    assign("current_values",init_out[['init_values']],envir = .GlobalEnv)
    assign("stored_vars",init_out[['stored_vars']],envir = .GlobalEnv)
    assign("hyperparameters",init_out[['hyperparameters']],envir = .GlobalEnv)
    assign("proposal_sigs",init_out[['proposal_sigs']],envir = .GlobalEnv)
    
    acceptance_rates=list()
    
    for(varname in varname_list)
    {
      acceptance_rates[[varname]]=rep(NA,M)
    }
    
    stored_parameters=list()
    for(varname in names(current_values))
    {
      stored_parameters[[varname]]=list()
      stored_parameters[[varname]][[1]]=current_values[[varname]]
    }
    
    stored_likelihoods=rep(NA,M)
    stored_likelihoods[1]=calculate_full_likelihood(stored_parameters,1)
    
    
    ###sample
    plot_fun(stored_parameters,1,mytitle='Initial Configuration',save_fig=save_fig,save_filename=paste(plot_dirname,'/initial_configuration.png',sep=''))
    start_index=2
  }
  
  pb <- progress_bar$new(
    format = " producing samples [:bar] :percent eta: :eta\n",
    total = M, clear = T, width= 80)
  print(data.frame('scale'=exp(current_values$logscale),'meantheta'=mean(current_values$theta),'meanbeta'=mean(current_values$beta)))
  
  for(jj in start_index:M)
  {
    current_seed=runif(1)*1e9 #for continuation purposes
    set.seed(current_seed)
    store=jj%%store_iter==0
    storej=floor(jj/store_iter)
    
    pb$tick()
    
    for(varname in varname_list)
    {
      out=update_vector(varname)
      
      global_update("current_values",varname,out$newvalue)
      
      if(store) {stored_parameters[[varname]][[storej]]=out$newvalue}
      
      acceptance_rates[[varname]][jj]=mean(out$accepts)
      
      if(latent_tf[[varname]])
      {
        normscale=sqrt(mean(rbind(current_values$z,current_values$w)^2))
        global_update("current_values","z",current_values$z/sqrt(mean(current_values$z^2)))
        global_update("current_values","w",current_values$w/sqrt(mean(current_values$w^2)))
        # global_update("current_values","mu",current_values$mu/normscale)
        # global_update("current_values","sigma",current_values$sigma/normscale)
      }
      
      if(update_sigma_tf[[varname]])
      {
        newsigma=update_sigma(varname)
        sig_varname=paste("sigma",varname,sep="_")
        if(store) {stored_parameters[[sig_varname]][[storej]]=newsigma}
        global_update("current_values",sig_varname,newsigma)
      }
    }
    
    new_is_spike=0 #update_is_spike()
    global_update("current_values",'is_spike',new_is_spike)
    if(store) {stored_parameters$is_spike[[storej]]=new_is_spike}
    
    
    if(ncluster>1)
    {
      assign("stored_vars",update_stored_vars(),envir = .GlobalEnv)
      
      newK=update_K()
      global_update("current_values","K_w",newK$K_w)
      global_update("current_values","K_z",newK$K_z)
      
      #store values for new Ks:
      global_update("current_values","mu",update_mu_clusters())
      global_update("current_values","sigma",update_sigma_clusters())
      global_update("current_values","lambda",update_lambda())
      global_update("current_values","omega",update_omega())
    }
    
    
    stored_likelihoods[jj]=calculate_full_likelihood(stored_parameters,jj)
    
    if(store) {      
      stored_parameters$K_w[[storej]]=current_values$K_w
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
      delta=ifelse(mean(acceptance_rates[[varname]][batch_size:jj],na.rm=T)>.44,min(.01,1/sqrt(jj)),-min(.01,1/sqrt(jj)))
      proposal_sigs[[varname]]=proposal_sigs[[varname]]*exp(delta)
    }
    assign("proposal_sigs",proposal_sigs,envir = .GlobalEnv)
    
    if(jj%%plot_iter==0)
    {
      # print('acceptance rates:')
      print(sapply(acceptance_rates,function(x) mean(x,na.rm=T)))
      print(data.frame('scale'=exp(current_values$logscale),'meantheta'=mean(current_values$theta),'meanbeta'=mean(current_values$beta)))
      plot_fun(stored_parameters,storej,mytitle=toString(jj),save_fig=save_fig,save_filename=paste(plot_dirname,'/iteration_',jj,'.png',sep=''))
      save(stored_parameters,stored_likelihoods,
           current_values,stored_vars,hyperparameters,
           proposal_sigs,acceptance_rates,current_seed,
           varname_list,update_sigma_tf,file=save_filename)
      
      myscale=c(exp(current_values$logscale))*(1-current_values$is_spike)
      wz_dist=myscale*euc_dist(current_values$z,current_values$w)
      bt_mat=outer(c(current_values$theta),c(current_values$beta),'+')
      preds=bt_mat-wz_dist
      print(preds)
      print(X)
      print(preds>0)
    }
  }
  return(stored_parameters)
}



M=1000
myseed=222
config_number=1
load_data=load_drv_data
ordinal=F
# load_data=load_spelling_data
# ordinal=FALSE
plot_iter=100
scaleterms=0
sampler_output=run_mcmc_sampler(1000,myseed,config_number,plot_iter=plot_iter,
                                load_data=load_data,ordinal=ordinal,save_fig=F,
                                ncluster=1,overwrite=T,return_output=T)


# est_scale=c()
# est_is_spike=c()
# for(scaleterm in c(0,1,2))
# {
#   print(paste('scaleterm:',scaleterm))
#   load_data=function() generate_simulated_data(scaleterm)
#   load_data()
#   sampler_output=run_mcmc_sampler(300,myseed,config_number,plot_iter=plot_iter,
#                                   load_data=load_data,ordinal=ordinal,save_fig=F,
#                                   ncluster=1,overwrite=T,return_output=T)
#   est_scale=c(est_scale,mean(exp(unlist(sampler_output$logscale))))
#   est_is_spike=c(est_is_spike,mean(unlist(sampler_output$is_spike)))
# }
# ordinal=F
# load_data=load_spelling_data
# ordinal=FALSE


# # 
# save_filename=file.path(paste('Saved_output/saved_output_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
# load(save_filename,verb=T)
# # rotated_parameters=procrustes_postprocess(stored_parameters,stored_likelihoods)
# # post_processed_parameters=cluster_relabeling(rotated_parameters)
# p=plot_latent_cluster(stored_parameters,10000,mytitle=toString(10000),save_fig=F,plot_pie=T)
# 
# 
# fivenum<-function(x)
# {
#   return(list('min'=min(x),'Q1'=quantile(x,.25),'median'=median(x),'Q3'=quantile(x,.75),'max'=max(x),'mean'=mean(x),'sd'=sd(x)))
# }
# average_response=apply(X,1,mean)
# 
# q0=data.frame(t(do.call('rbind',lapply(1:ncluster,function(clusti) apply(X[current_values$K_z==clusti,],2,mean)))))
# q1=data.frame(t(do.call('rbind',lapply(1:ncluster,function(clusti) fivenum(current_values$theta[current_values$K_z==clusti])))))
# q2=data.frame(t(do.call('rbind',lapply(1:ncluster,function(clusti) fivenum(average_response[current_values$K_z==clusti])))))
# names(q1)=c('cluster1','cluster2','cluster3')
# names(q2)=c('cluster1','cluster2','cluster3')
# xtable(q0)
# xtable(q1)
# xtable(q2)

# save_image(p,'../Reports/Figures/simulated_cluster.png')

# source('postprocessing_code.R')
# rotated_parameters=procrustes_postprocess(stored_parameters,stored_likelihoods)
# post_processed_parameters=cluster_relabeling(rotated_parameters)
# p=plot_latent_cluster(post_processed_parameters,10000,mytitle=toString(10000),save_fig=F,plot_pie=T)
# print(p)
# save_image(p,'../Reports/Figures/abortion_cluster3.png')
# 
# 
# qclust=data.frame(t(do.call('rbind',lapply(1:ncluster,function(clusti) apply(X[current_values$K_z==clusti,],2,mean)))))
# names(qclust)=c('cluster1','cluster2','cluster3')
# qclust=rbind(qclust,'average'=apply(qclust,2,mean))
# xtable(qclust)

# print(sapply(acceptance_rates,function(x) mean(x,na.rm=T)))


# # 
# save_filename=file.path(paste('Saved_output/saved_output_config_',config_number,'_seed_',myseed,'_data_',dataname,sep=''))
# load(save_filename,verb=T)
# # rotated_parameters=procrustes_postprocess(stored_parameters,stored_likelihoods)
# # post_processed_parameters=cluster_relabeling(rotated_parameters)
# p=plot_latent_cluster(stored_parameters,10000,mytitle=toString(10000),save_fig=F,plot_pie=T)
# 
# 
# fivenum<-function(x)
# {
#   return(list('min'=min(x),'Q1'=quantile(x,.25),'median'=median(x),'Q3'=quantile(x,.75),'max'=max(x),'mean'=mean(x),'sd'=sd(x)))
# }
# average_response=apply(X,1,mean)
# 
# q0=data.frame(t(do.call('rbind',lapply(1:ncluster,function(clusti) apply(X[current_values$K_z==clusti,],2,mean)))))
# q1=data.frame(t(do.call('rbind',lapply(1:ncluster,function(clusti) fivenum(current_values$theta[current_values$K_z==clusti])))))
# q2=data.frame(t(do.call('rbind',lapply(1:ncluster,function(clusti) fivenum(average_response[current_values$K_z==clusti])))))
# names(q1)=c('cluster1','cluster2','cluster3')
# names(q2)=c('cluster1','cluster2','cluster3')
# xtable(q0)
# xtable(q1)
# xtable(q2)

# save_image(p,'../Reports/Figures/simulated_cluster.png')

# source('postprocessing_code.R')
# rotated_parameters=procrustes_postprocess(stored_parameters,stored_likelihoods)
# post_processed_parameters=cluster_relabeling(rotated_parameters)
# p=plot_latent_cluster(post_processed_parameters,10000,mytitle=toString(10000),save_fig=F,plot_pie=T)
# print(p)
# save_image(p,'../Reports/Figures/abortion_cluster3.png')
# 
# 
# qclust=data.frame(t(do.call('rbind',lapply(1:ncluster,function(clusti) apply(X[current_values$K_z==clusti,],2,mean)))))
# names(qclust)=c('cluster1','cluster2','cluster3')
# qclust=rbind(qclust,'average'=apply(qclust,2,mean))
# xtable(qclust)

# print(sapply(acceptance_rates,function(x) mean(x,na.rm=T)))
