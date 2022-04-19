######################################################
##Background: in the script, we firstly take one site from ETH(Chamau, GRA, switerland)
##to do the data tidy for the further analysis
Sys.setenv(tz='UTC')
##function
##load the library
library(ncdf4)
library(lubridate)
library(phenopix)
library(plyr)
##
FLUXNETcodes<-c('CH-Cha',"AT-Neu",
                "DE-Hai")
#"US-Var","US-Wkg","US-Ha1","CA-Oas"
FLUXpath<-'M:/data/DataStructureMDI/DATA/Incoming/Fluxnet/berkeley_012016/Data/HH/2016_11/'
n<-length(FLUXNETcodes)
for(i in 3:n){
  fname<-list.files(path = FLUXpath, pattern=FLUXNETcodes[i],
                    full.names = TRUE)
  fname_short<-paste0(FLUXNETcodes[i],' ','Fluxes&Meteos')
  if (length(fname) > 0){
    print(paste('###################### Analyzing ', fname_short))
    
    nc <- nc_open(fname)
    
    ##extract the variables from the nc file
    #Date and Time
    Year <- ncvar_get(nc, "year")
    Month<-ncvar_get(nc,'month')
    Day<-ncvar_get(nc,'day')
    Hour<-ncvar_get(nc, "hour")
    Min<-ncvar_get(nc, "minute")
    Night<-ncvar_get(nc, "NIGHT")
    
    ###Carbon Fluxes
    GPP_NT_VUT_MEAN<- ncvar_get(nc, "GPP_NT_VUT_MEAN")
    GPP_NT_VUT_SE<- ncvar_get(nc, "GPP_NT_VUT_SE")
    NEE_VUT_MEAN<- ncvar_get(nc,"NEE_VUT_MEAN")
    NEE_VUT_SE<- ncvar_get(nc,"NEE_VUT_SE")
    NEE_VUT_MEAN_QC<- ncvar_get(nc,"NEE_VUT_MEAN_QC")
    RECO_NT_VUT_MEAN<- ncvar_get(nc,"RECO_NT_VUT_MEAN")
    RECO_NT_VUT_SE<- ncvar_get(nc,"RECO_NT_VUT_SE")
    
    ##Engery Fluxex
    G_F_MDS<-ncvar_get(nc,"G_F_MDS")
    G_F_MDS_QC<-ncvar_get(nc,"G_F_MDS_QC")
    LE_F_MDS<-ncvar_get(nc,"LE_F_MDS")
    LE_F_MDS_QC<-ncvar_get(nc,"LE_F_MDS_QC")
    H_F_MDS<-ncvar_get(nc,"H_F_MDS")
    H_F_MDS_QC<-ncvar_get(nc,"H_F_MDS_QC")
    
    ##Meteos
    SW_IN_F_MDS<-ncvar_get(nc,"SW_IN_F_MDS")
    SW_IN_F_MDS_QC<-ncvar_get(nc,"SW_IN_F_MDS_QC")
    SW_OUT<-ncvar_get(nc,"SW_OUT")
    SW_OUT_QC<-ncvar_get(nc,"SW_OUT_QC")
    PPFD_IN<-ncvar_get(nc,"PPFD_IN")
    PPFD_IN_QC<-ncvar_get(nc,"PPFD_IN_QC")
    PPFD_OUT<-ncvar_get(nc,"PPFD_OUT")
    PPFD_OUT_QC<-ncvar_get(nc,"PPFD_OUT_QC")
    TA_F_MDS<-ncvar_get(nc,"TA_F_MDS")
    TA_F_MDS_QC<-ncvar_get(nc,"TA_F_MDS")
    VPD_F_MDS<-ncvar_get(nc,"VPD_F_MDS")
    VPD_F_MDS_QC<-ncvar_get(nc,"VPD_F_MDS_QC")
    P_F<-ncvar_get(nc,"P_F")
    P_F_QC<-ncvar_get(nc,"P_F_QC")
    SWC_F_MDS_1<-ncvar_get(nc,"SWC_F_MDS_1")
    SWC_F_MDS_1_QC<-ncvar_get(nc,"SWC_F_MDS_1_QC")
    SWC_F_MDS_2<-ncvar_get(nc,"SWC_F_MDS_2")
    SWC_F_MDS_2_QC<-ncvar_get(nc,"SWC_F_MDS_2_QC")
    SWC_F_MDS_3<-ncvar_get(nc,"SWC_F_MDS_3")
    SWC_F_MDS_3_QC<-ncvar_get(nc,"SWC_F_MDS_3_QC")
    #ncatt_get(nc,"G_F_MDS")
    ##time setting
    timeStampChar <-paste0(Year,'-',Month,'-',Day,' ',Hour,':',Min)
    Date<-ymd_hm(timeStampChar)
    
    #produce the new datasets
    Data.F<-data.frame(Year=Year,Month=Month,Day=Day,Hour=Hour,Date=Date,Night=Night,
                       GPP=GPP_NT_VUT_MEAN,GPP_SE=GPP_NT_VUT_SE,
                       NEE=NEE_VUT_MEAN,NEE_SE=NEE_VUT_SE,NEE_QC=NEE_VUT_MEAN_QC,
                       RECO=RECO_NT_VUT_MEAN,RECO_SE=RECO_NT_VUT_SE,
                       G=G_F_MDS,G_QC=G_F_MDS_QC,
                       LE=LE_F_MDS,LE_QC=LE_F_MDS_QC,
                       H=H_F_MDS,H_QC=H_F_MDS_QC,
                       SW_IN=SW_IN_F_MDS,SW_IN_QC=SW_IN_F_MDS_QC,
                       SW_OUT=SW_OUT,SW_OUT_QC=SW_OUT_QC,
                       PPFD_IN=PPFD_IN,PPFD_IN_QC=PPFD_IN_QC,
                       PPFD_OUT=PPFD_OUT,PPFD_OUT_QC=PPFD_OUT_QC,
                       TA=TA_F_MDS,TA_QC=TA_F_MDS_QC,
                       VPD=VPD_F_MDS,VPD_QC=VPD_F_MDS_QC,
                       P=P_F,P_QC=P_F_QC,
                       SWC_1=SWC_F_MDS_1,SWC_1_QC=SWC_F_MDS_1_QC,
                       SWC_2=SWC_F_MDS_2,SWC_2_QC=SWC_F_MDS_2_QC,
                       SWC_3=SWC_F_MDS_3,SWC_3_QC=SWC_F_MDS_3_QC
                       )
    ##Calculate the daily sum of Precipitation(P)
    Data.F$rDate<-as.POSIXct(strptime(Data.F$Date,format = "%Y-%m-%d"))
    P_d<-ddply(.data = Data.F,.(rDate),summarise,P_d=sum(P))
    ##Take out the data that in the daytime-->Night==0
    Data.F<-Data.F[(Data.F$Night==0),]
    ##->Selected high quality data for NEE-->NEE_QC==0
    Data.F<-Data.F[Data.F$NEE_QC==0|Data.F$NEE_QC==1,]
    
    ##1. Average to Daily
    Data.F$NEP<-Data.F$NEE*(-1)
    ##note:SWC was averaged using first three layers data
    Data.F$SWC_Mean<-(Data.F$SWC_1+Data.F$SWC_2+Data.F$SWC_3)/3
    Data.F_daily<-ddply(.data = Data.F,.(rDate),summarise,
                        GPP_d=mean(GPP,na.rm = T),NEP_d=mean(NEP,na.rm = T),RECO_d=mean(RECO,na.rm = T),
                        G_d=mean(G,na.rm = T),LE_d=mean(LE,na.rm = T),H_d=mean(H,na.rm = T),
                        SW_IN_d=mean(SW_IN,na.rm = T),SW_OUT_d=mean(SW_OUT,na.rm = T),
                        PPFD_IN_d=mean(PPFD_IN,na.rm = T),PPFD_OUT_d=mean(PPFD_OUT,na.rm = T),
                        TA_d=mean(TA,na.rm = T),VPD_d=mean(VPD,na.rm = T),
                        SWC_d=mean(SWC_Mean,na.rm = T))
    ##a.add NAs to make full time series
    #find the minimum time and max time
    minT<-min(Data.F_daily$rDate,na.rm=T)
    maxT<-max(Data.F_daily$rDate,na.rm=T)
    #create full time series
    tm=60*60*24
    tmin<-as.POSIXct(strptime(minT,format = '%Y-%m-%d'));tmax<-as.POSIXct(strptime(maxT,format = '%Y-%m-%d'))
    TS<-seq(tmin, tmax, by = "day")
    TS<-data.frame(rDate=as.POSIXct(TS))
    #fill gaps in time series with NA 
    temp<-merge(TS,Data.F_daily,by.x='rDate',all.x=T)
    #also add the sum of daily Preci into the datasets
    P_d<-merge(TS,P_d,by.x = "rDate")
    ##produced the full time series of daily data
    Data.F_daily<-cbind(temp,P_d=P_d$P_d)
    
    #b.remove the NAs at the start of time series (keep the time series only when they including TA, SW_In and P)
    remove_headNAs<-function(df,index_name){
      # df<-Data.F_daily
      # index_name<-'SW_IN_d'
      
      N=length(df[,index_name])
      pos_1NoNA<-min(which(!is.na(df[,index_name])))
      df<-df[pos_1NoNA:N,]
      return(df)
    }
    Data.F_daily<-remove_headNAs(Data.F_daily,'TA_d')
    Data.F_daily<-remove_headNAs(Data.F_daily,'SW_IN_d')
    Data.F_daily<-remove_headNAs(Data.F_daily,'P_d')
    
    #c.remove the outliers using spline methods(exclude the data>2*sd)
    source(file=paste0("M:/people/yluo/Ecosystem_stability/code/Sites_based_analysis/test/functions/SplineFilter_forFluxnets.R"))
    #only choose most important variables
    #->  "rDate"      "GPP_d"      "NEP_d"      "RECO_d"     "G_d"        "LE_d"       "H_d"       
    # "SW_IN_d"        "TA_d"       "VPD_d"      "SWC_d"     "P_d" 
    Sel_Vars<-c("GPP_d","NEP_d","RECO_d","G_d","LE_d","H_d","SW_IN_d","TA_d","VPD_d","SWC_d")
    Data.F_daily_filter<-apply(Data.F_daily[,Sel_Vars],2,SplineFilter)
    Data.F_daily_filter<-data.frame(rDate=Data.F_daily$rDate,Data.F_daily_filter)
    #Precipitation is random variables(not continuous), should not use the filter to filter
    Data.F_daily_filter$P_d<-Data.F_daily$P_d
    
    #save the tidyed data  
    main_path<-'M:/people/yluo/Ecosystem_stability/Data_for_use/Daily_filter_data/'
    outRdatapath<-paste0(main_path, FLUXNETcodes[i], "_Daily.RDA")

    print(paste("############################ Saving ", fname_short,'which has been filtered with spline'))
    save(Data.F_daily_filter, file = outRdatapath) 
  }
}
