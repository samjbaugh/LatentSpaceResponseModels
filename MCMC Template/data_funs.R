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
  raw_data=read.table('../Data/spelling.dat')
  assign("X",VerbAgg,envir=.GlobalEnv)
  assign('nz',dim(X)[1],envir=.GlobalEnv)
  assign('nw',dim(X)[2],envir=.GlobalEnv)
}