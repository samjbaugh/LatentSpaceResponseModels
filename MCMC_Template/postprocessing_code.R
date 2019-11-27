require('expm')

procrustes_postprocess<-function(stored_parameters,stored_likelihoods)
{
  nz=dim(stored_parameters$z[[1]])[1]
  nw=dim(stored_parameters$w[[1]])[1]
  M=sum(!is.na(stored_likelihoods))
  stored_likelihoods=stored_likelihoods[1:M]
  imax=which(stored_likelihoods==max(stored_likelihoods,na.rm=T))
  wzmax=rbind(stored_parameters$z[[imax]],stored_parameters$w[[imax]])
  mean_pos=apply(wzmax,2,mean)
  wz0=sweep(wzmax,2,mean_pos)
  
  matched_parameters=stored_parameters
  
  for(ii in 1:M)
  {
    #center first
    wz=rbind(stored_parameters$z[[ii]],stored_parameters$w[[ii]])
    mean_pos=apply(wz,2,mean)
    wz_centered=sweep(wz,2,mean_pos)
    mu_centered=sweep(stored_parameters$mu[[ii]],2,mean_pos)
    
    rot_mat=t(wz0)%*%wz_centered%*%solve(sqrtm(t(wz_centered)%*%wz0%*%t(wz0)%*%wz_centered))
    wz_matched=wz_centered%*%rot_mat
    mu_matched=mu_centered%*%rot_mat
    
    matched_parameters$z[[ii]]=wz_matched[1:nz,]
    matched_parameters$w[[ii]]=wz_matched[(nz+1):(nz+nw),]
    matched_parameters$mu[[ii]]=mu_matched
  }
  return(matched_parameters)
}

cluster_relabeling<-function(parameterlist,order_var='sigma')
{
  ncluster=length(unique(parameterlist$K_z[[1]]))
  relabeled_parameters=parameterlist
  for(ii in 1:length(parameterlist))
  {
    order_stat=relabeled_parameters[[order_var]][[ii]]
    new_ordering=order(order_stat)
    relabeled_parameters$sigma[[ii]]=relabeled_parameters$sigma[[ii]][new_ordering]
    relabeled_parameters$mu[[ii]]=relabeled_parameters$mu[[ii]][new_ordering,]
    relabeled_parameters$K_w[[ii]]=factor(relabeled_parameters$K_w[[ii]],levels=1:ncluster,labels=new_ordering)
    relabeled_parameters$K_z[[ii]]=factor(relabeled_parameters$K_z[[ii]],levels=1:ncluster,labels=new_ordering)
  }
  
  return(relabeled_parameters)
}
# # 
# # load_charity_data()
# # load(file='./Saved_output/saved_output_config_1_seed_555_data_spelling',verb=T)
# # matched_parameters=procrustes_postprocess(stored_parameters,stored_likelihoods)
# # init_out=initialize_sampler(2,ordinal=F)
# 
# plot_dirname=file.path(paste('Images/plots_config_',1,'_seed_',555,'_data_','spelling',sep=''))
# plot_fun=plot_latent_cluster
# for(jj in seq(100,2000,by=100))
# {
#   plot_fun(matched_parameters,jj,mytitle=toString(jj),save_fig=T,save_filename=paste(plot_dirname,'/iteration_',jj,'_matched.png',sep=''))
# }
# 
# plot_fun(matched_parameters,1565,mytitle=toString(1565),save_fig=T,save_filename=paste(plot_dirname,'/iteration_',1565,'_matched.png',sep=''))
# plot_fun(stored_parameters,1565,mytitle=toString(1565),save_fig=T,save_filename=paste(plot_dirname,'/iteration_',1565,'.png',sep=''))
