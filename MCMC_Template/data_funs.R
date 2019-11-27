generate_simulated_data<-function(scale_term)
{
  set.seed(123)
  nw=4
  nz=10
  theta=rnorm(nz)
  beta=rnorm(nw)
  ztrue=matrix(rnorm(nz*2),nz,2)
  ztrue=ztrue/c(sqrt(mean(ztrue^2)))
  wtrue=matrix(rnorm(nw*2),nw,2)
  wtrue=wtrue/c(sqrt(mean(wtrue^2)))
  bt_mat=outer(theta,beta,'+')
  wz_dist=euc_dist(ztrue,wtrue)
  scale_term=scale_term
  X_probs=sigmoid(bt_mat-scale_term*wz_dist)
  X=matrix(rbinom(n=prod(dim(X_probs)),size=1,prob=X_probs),dim(X_probs))
  assign("X",X,envir=.GlobalEnv)
  assign("true_theta",theta,envir=.GlobalEnv)
  assign("true_beta",beta,envir=.GlobalEnv)
  assign("dataname","simulated_large",envir=.GlobalEnv)
  assign('nz',dim(X)[1],envir=.GlobalEnv)
  assign('nw',dim(X)[2],envir=.GlobalEnv)
}

load_abortion_data<-function()
{
  raw_data=read.table('../Data/abortion.txt')
  
  assign("X",raw_data,envir=.GlobalEnv)
  assign("dataname","abortion",envir=.GlobalEnv)
  assign('nz',dim(X)[1],envir=.GlobalEnv)
  assign('nw',dim(X)[2],envir=.GlobalEnv)
}

load_spelling_data<-function()
{
  raw_data=read.table('../Data/spelling.dat')
  
  assign("X",raw_data[,2:5],envir=.GlobalEnv)
  assign("dataname","spelling",envir=.GlobalEnv)
  assign("gender",factor(raw_data[,1],levels=c(0,1),labels=c('female','male')),envir=.GlobalEnv)
  assign('nz',dim(X)[1],envir=.GlobalEnv)
  assign('nw',dim(X)[2],envir=.GlobalEnv)
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

load_big5_data<-function()
{
  rawdata=read.table('../Data/big5.csv', header = TRUE, fill = TRUE)
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

load_verbal_agression_data<-function()
{
  require('irtrees')
  data("VerbAgg3")
  
  assign("X",VerbAgg3[,3:dim(VerbAgg3)[2]]-1,envir=.GlobalEnv)
  assign('nz',dim(X)[1],envir=.GlobalEnv)
  assign('nw',dim(X)[2],envir=.GlobalEnv)
  
  assign("dataname","verbagg",envir=.GlobalEnv)
  
  ordinal_categories=unique(c(as.matrix(X)))
  assign("ordinals",ordinal_categories[!is.na(ordinal_categories)],envir=.GlobalEnv)
  assign('minord',min(ordinals),envir=.GlobalEnv)
  assign('maxord',max(ordinals),envir=.GlobalEnv)
  assign('ntau',length(ordinals),envir=.GlobalEnv)
}