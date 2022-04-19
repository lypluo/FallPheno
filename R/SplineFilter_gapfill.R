##########################
#function used to filter the data
##########################
library(zoo)
library(phenopix)
SplineFilter_thengap_f<-function(ts_ori){
  # ts_ori<-df.proc.VIs[,2]
  ##do primary gapfilling in order to enable SplineFit function afterwards
  ##do not use SSA method as sometimes it extrapolates too much
  # ts_ori<-df.proc.VIs$EVI
  #for Vegetation index, set lower limit <0
  hard_lower_limit<-0
  #hard removing
  ts_ori[ts_ori<hard_lower_limit]<-NA
  #remove based on the mean and sd of the times series(based on 2*sd):
  mean_ts<-mean(ts_ori,na.rm=T)
  sd_ts<-sd(ts_ori,na.rm = T)
  ts_ori[ts_ori<c(mean_ts-2*sd_ts)]<-NA
  pos_NA_ori<-which(is.na(ts_ori))
  #in case some NA at the beginning of the time series,use na.fill for the beginning of gapfilling
  temp1<-na.fill(ts_ori[1:150],fill='extend') ##here temporily set the 150 day as the maximum flling gap
  temp2<-na.spline(ts_ori[151:length(ts_ori)])
  ts_new<-c(temp1,temp2)
  
  ##########transfer the date format to zoo
  pos_noise<-c()
  for(i in 1:20){
    ts_new<-zoo(ts_new,order.by = index(1:length(ts_new)))
    ####delete the data which bigger than (SplineFit-sd*2,SplineFit+sd*4)
    fitResult<-SplineFit(ts_new,nrep=100,df.factor = 0.05)
    residuals<-as.numeric(ts_new)-as.vector(fitResult$fit$predicted)
    sd.res<-sd(residuals,na.rm=TRUE)
    selectcriterion_above<-mean(residuals)+2*sd.res
    selectcriterion_below<-mean(residuals)-2*sd.res
    noise<-residuals[residuals>=selectcriterion_above|residuals<=selectcriterion_below]
    matchNum<-match(noise,residuals)
    pos_noise<-c(pos_noise,matchNum)
    if(length(matchNum)==0){
      break
    }
    ts_new[matchNum]<-NA
    ts_NA1<-as.numeric(ts_new)
    ###gapfilling with na.fill function
    temp1<-na.fill(ts_NA1,fill='extend')
    ts_new<-temp1
  }
  #set all the outliers(NA_ori) to NA
  ts_new[pos_NA_ori]<-NA
  return(ts_new)
}

