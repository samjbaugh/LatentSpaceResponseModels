load_charity_data<-function()
{
  rawdata=read.table('../Data/charity.dat')
  X=rawdata[!apply(is.na(rawdata),1,any),]
  return(list('data'=X,'dataname'='charity'))
}