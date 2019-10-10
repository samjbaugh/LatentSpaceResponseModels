load_spelling_data<-function()
{
  raw_data=read.table('../Data/spelling.dat')
  return(list('data'=raw_data[,2:5],'dataname'='spelling','gender'=factor(raw_data[,1],levels=c(0,1),labels=c('female','male'))))
}
