load_spelling_data<-function()
{
  raw_data=read.table('../Data/spelling.dat')
  
  assign("X",raw_data[,2:5],envir=.GlobalEnv)
  assign("dataname","spelling",envir=.GlobalEnv)
  assign("gender",factor(raw_data[,1],levels=c(0,1),labels=c('female','male')),envir=.GlobalEnv)
}

load_charity_data<-function()
{
  rawdata=read.table('../Data/charity.dat')
  assign("X",rawdata[!apply(is.na(rawdata),1,any),],envir=.GlobalEnv)
  assign('nz',dim(X)[1],envir=.GlobalEnv)
  assign('nw',dim(X)[2],envir=.GlobalEnv)
  
  assign("dataname","charity",envir=.GlobalEnv)
  
  ordinal_categories=unique(c(as.matrix(X)))
  assign("ordinals",ordinal_categories[!is.na(ordinal_categories)],envir=.GlobalEnv)
  assign('minord',min(ordinals),envir=.GlobalEnv)
  assign('maxord',max(ordinals),envir=.GlobalEnv)
  assign('ntau',length(ordinals),envir=.GlobalEnv)
}

simulate_data_longitudinal<-function(myseed)
{
  assign("nz",100,envir=.GlobalEnv)
  assign("nw",4,envir=.GlobalEnv)
  
  set.seed(myseed)
  true_sigma_z1=1
  true_sigma_z2=1
  true_sigma_w=1
  # sd_init=sqrt(100)
  true_theta1=matrix(rnorm(nz,0,1),nz,1)
  true_theta2=true_theta1+.5 #rnorm(nz,0,1)
  true_beta=matrix(rnorm(nw,0,1),nw,1)
  true_z2=matrix(cbind(rnorm(nz,2,true_sigma_z1),rnorm(nz,2,true_sigma_z1)),nz,2)
  true_z1=matrix(cbind(rnorm(nz,-2,true_sigma_z1),rnorm(nz,-2,true_sigma_z1)),nz,2)
  colnames(true_z1)<-c('coord1','coord2')
  colnames(true_z2)<-c('coord1','coord2')
  true_w=matrix(cbind(rnorm(nw,2*c(-2,-2,2,2),true_sigma_w),rnorm(nw,2*c(-2,-2,2,2),true_sigma_w)),nw,2)
  colnames(true_w)<-c('coord1','coord2')
  
  bt_mat1=outer(c(true_theta1),c(true_beta),'+')
  bt_mat2=outer(c(true_theta2),c(true_beta),'+')
  
  wz_dist1=euc_dist(true_z1,true_w)
  wz_dist2=euc_dist(true_z2,true_w)
  
  X1_probs=sigmoid(-bt_mat1-wz_dist1)
  X2_probs=sigmoid(-bt_mat2-wz_dist2)
  # set.seed(123)
  
  X1=matrix(rbinom(n=prod(dim(X1_probs)),size=1,prob=X1_probs),dim(X1_probs))
  X2=matrix(rbinom(n=prod(dim(X2_probs)),size=1,prob=X2_probs),dim(X2_probs))
  plot_simulated=T
  if(plot_simulated)
  {
    #for plotting simulated data
    simulated=list()
    simulated$z1[[1]]=true_z1
    simulated$z2[[1]]=true_z2
    simulated$w[[1]]=true_w
    
    print(plot_latent_longitudinal(simulated,1,mytitle="Simulated Data True Latent Space",cluster=F)) 
  }
  
  assign("X1",X1,envir=.GlobalEnv)
  assign("X2",X2,envir=.GlobalEnv)
  assign("dataname","longitudinalsim",envir=.GlobalEnv)
}

load_big5_data<-function()
{
  rawdata=read.table('./Data/big5.csv', header = TRUE, fill = TRUE)
  rawdata <- select(rawdata, E1:E10) #only extraversion questions
  rawdata[rawdata == 0] <- NA
  rawdata <- rawdata[!apply(is.na(rawdata),1,any),]
  
  #only taking a sample
  samp <- sample(1:nrow(rawdata), 50)
  rawdata <- rawdata[samp,]

  assign("X", rawdata ,envir=.GlobalEnv)
  assign('nz',dim(X)[1],envir=.GlobalEnv)
  assign('nw',dim(X)[2],envir=.GlobalEnv)
  
  assign("dataname","BIG5",envir=.GlobalEnv)
  
  ordinal_categories=unique(c(as.matrix(X)))
  assign("ordinals",ordinal_categories[!is.na(ordinal_categories)],envir=.GlobalEnv)
  assign('minord',min(ordinals),envir=.GlobalEnv)
  assign('maxord',max(ordinals),envir=.GlobalEnv)
  assign('ntau',length(ordinals),envir=.GlobalEnv)
}