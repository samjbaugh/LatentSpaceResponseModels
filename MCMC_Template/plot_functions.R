
plot_latent=function(stored_parameters,store_index,mytitle="",save_fig=F,save_filename="")
{
  z=data.frame(stored_parameters$z[[store_index]]) #data.frame(matrix(stored_parameters$z[M,],nz,2))
  names(z)<-c('coord1','coord2')
  
  w=data.frame(stored_parameters$w[[store_index]]) #data.frame(matrix(stored_parameters$w[M,],nw,2))
  names(w)<-c('coord1','coord2')
  w$name=sapply(1:nw,toString)
  
  p0<-ggplot()+geom_point(aes(x=coord1,y=coord2,col='z'),z,cex=2)+
    xlab('coordinate 1')+ylab('coordinate 2')+
    ggtitle(paste('Latent space sample at k=',mytitle,sep=''))+
    geom_text(aes(x=coord1,y=coord2,label=name,col='w'),hjust=2,vjust=2,data=w,fontface="bold")+
    coord_fixed()
  
  if(save_fig)
  {
    png(save_filename)
    print(p0)
    dev.off()
  }else
  {
    print(p0)
  }
}


plot_latent_cluster <- function(stored_parameters,store_index,mytitle="",save_fig=F,plot_pie=F,save_filename="")
{
  z=data.frame(stored_parameters$z[[store_index]]) #data.frame(matrix(stored_parameters$z[M,],nz,2))
  names(z)<-c('coord1','coord2')
  
  w=data.frame(stored_parameters$w[[store_index]]) #data.frame(matrix(stored_parameters$w[M,],nw,2))
  names(w)<-c('coord1','coord2')
  
  temp_z=sapply(stored_parameters$K_z[[store_index]],function(x) paste(toString(x)))
  z$K=as.factor(temp_z)
  
  temp_w=sapply(stored_parameters$K_w[[store_index]],function(x) paste(toString(x)))
  w$K=as.factor(temp_w)
  w$name=sapply(1:nw,toString)
  
  if(plot_pie)
  {
    total=length(stored_parameters$K_w)
    clust_dist <- function(wnum) {return(sapply(1:ncluster,function(iclust) mean(sapply(stored_parameters$K_w[burn_in:total],function(x) x[wnum]==iclust))))}
    clust_totals <- sapply(1:nw,clust_dist)
    for(iclust in 1:ncluster)
    {
      w[[paste('cluster',iclust,sep='')]]=clust_totals[iclust,]
    }
  }
  
  mu=data.frame(stored_parameters$mu[[store_index]])
  colnames(mu)<-c('coord1','coord2')
  mu$K=sapply(1:ncluster,function(x) paste(toString(x)))
  mu$r=stored_parameters$sigma[[store_index]]
  
  p0<-ggplot()+geom_point(aes(x=coord1,y=coord2,col=K),z,cex=2)+
    xlab('coordinate 1')+ylab('coordinate 2')+
    ggtitle(paste('Latent space sample at k=',mytitle,sep=''))+
    geom_point(aes(x=coord1,y=coord2),pch=8,col='black',data=mu)+
    geom_circle(aes(x0=coord1,y0=coord2,r=r),data=mu)+
    geom_text(aes(x=coord1,y=coord2,label=name),hjust=3,vjust=2,data=w,fontface="bold")+
    coord_fixed()
  if(plot_pie){
  p0<-p0+geom_scatterpie(aes(x=coord1,y=coord2,r=.08),data=w,cols=sapply(1:ncluster,function(iclust) paste('cluster',iclust,sep='')))}
  p0
  if(save_fig)
  {
    png(save_filename)
    print(p0)
    dev.off()
  }else
  {
    print(p0)
  }
}

plot_latent_ordinal_cluster <- function(stored_parameters,store_index,mytitle="",save_fig=F,plot_pie=F,save_filename="")
{
  z=data.frame(stored_parameters$z[[store_index]]) #data.frame(matrix(stored_parameters$z[M,],nz,2))
  w=data.frame(stored_parameters$w[[store_index]]) #data.frame(matrix(stored_parameters$w[M,],nw,2))
  # nw=dim(w)[1]
  # nz=dim(z)[1]
  names(z)<-c('coord1','coord2')
  names(w)<-c('coord1','coord2')
  
  temp_z=sapply(stored_parameters$K_z[[store_index]],function(x) paste(toString(x)))
  z$K=as.factor(temp_z)
  
  temp_w=sapply(stored_parameters$K_w[[store_index]],function(x) paste(toString(x)))
  w$K=as.factor(temp_w)
  w$wname=c(rep('no',nw),rep('perhaps',nw),rep('yes',nw))#c(sapply(1:ntau,function(ii) rep(paste("w_k=",ii,sep=''),nw)))
  
  mu=data.frame(stored_parameters$mu[[store_index]])
  colnames(mu)<-c('coord1','coord2')
  mu$wz="mu"
  mu$K=sapply(1:ncluster,function(x) paste(toString(x)))
  mu$r=stored_parameters$sigma[[store_index]]
  
  p0<-ggplot()+
    geom_point(aes(x=coord1,y=coord2,col=K,pch='z'),data=z,cex=2)+
    geom_point(aes(x=coord1,y=coord2,col=K,pch=wname),data=w,cex=2)+
    geom_point(aes(x=coord1,y=coord2,col=K,pch='mu'),data=mu)+
    xlab('coordinate 1')+ylab('coordinate 2')+
    ggtitle(paste('Latent space sample at M=',mytitle,sep=''))+
    geom_circle(aes(x0=coord1,y0=coord2,col=K,r=r),data=mu)+
    scale_shape_manual(values=c("z"=16,"no"=0,"perhaps"=11,"yes"=8,"w_k=4"=11,
                                     "w_k=5" = 12, "mu"=4))
  # p1=p0+geom_point(aes(x=current_values[['mu_z']][,1],y=current_values[['mu_z']][,2],col=c("1_z","2_z","3_z","4_z","5_z"),pch="cluster mu z"),cex=4)
  if(save_fig)
  {
    png(save_filename)
    print(p0)
    dev.off()
  }else
  {
    print(p0)
  }
}

save_image<-function(myimage,imfilename)
{
  png(imfilename, width=3200, height=3200, res=600)
  print(myimage)
  dev.off()
  system(paste("convert -trim",imfilename,imfilename))
}

