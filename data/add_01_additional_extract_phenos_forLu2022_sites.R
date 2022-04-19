#-------------------------------
#Aim：To extract the phenos from the VIs from Xinchen's used site in Lu et al., 2022
#To compare if my extracted VIs has significant difference compared Lu et al., 2022
#-------------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
#--------------------------------------
#1)check the Lu's analysis sites:
#--------------------------------------
#LU2022 original sites
load(file = paste0("./data/data_prep/Add_analysis_based_Lu2022/Data_from_Lu2022.RDA"))
t1<-unique(df_merge$sitename)
#sites in Lu2022 we can access for VIs from Walther et al., 2021 datasets
load(paste0("./data/data_prep/Add_analysis_based_Lu2022/MODIS_VIs_fromWalther2021_forLU2022_sites.rda"))
t2<-unique(df.VIs$sitename)
#the sites are not aviable：
setdiff(t1,t2)

#------------------------------------
#2) start to extract the Phenos from VIs
#-----------------------------------
#-----2a) first check time series of data 
view_ts<-function(df,VI_name){
  # df<-df.VIs
  # VI_name<-"EVI"
  
  site.names<-unique(df$sitename)
  plot.all_doy<-c()
  plot.all_date<-c()
  for(i in 1: length(site.names)){
    df.proc<- df %>%
      mutate(Year=year(Date))%>%
      filter(sitename==site.names[i])
    #
    df.temp<-df.proc[,c("sitename","Year","Date",VI_name)]
    names(df.temp)<-c("sitename","Year","Date","VI")
    df.temp<-df.temp %>% mutate(doy=lubridate::yday(Date))
    #
    df.temp$Year<-as.factor(df.temp$Year)
    #plotting
    y.range<-range(df.temp$VI,na.rm=T)
    p_plot1<-ggplot(data = df.temp,aes(x=doy,y=VI,col=Year))+
      # geom_line()+
      geom_point()+
      ylab(VI_name)+
      annotate(geom = "text",x=250,y.range[1],
               label=paste0(df.temp$PFT,"-->",site.names[i]))
    plot.all_doy[[i]]<-p_plot1
    
    p_plot2<-ggplot(data = df.temp,aes(x=Date,y=VI,col=Year))+
      # geom_line()+
      geom_point()+
      ylab(VI_name)+
      annotate(geom = "text",x=as.Date("2010-06-30"),y.range[1],
               label=paste0(df.temp$PFT,"-->",site.names[i]))
    plot.all_date[[i]]<-p_plot2
  }
  #merge
  plot.all<-list(plot.all_date=plot.all_date,
                 plot.all_doy=plot.all_doy)
  return(plot.all)
}

#for VI:EVI
EVI.plots<-view_ts(df.VIs,"EVI")
#plotting
#Date
EVI.plots$plot.all_date
#DOY
# EVI.plots$plot.all_doy

#-----2b).Remove the outliers of ts---------
source("./R/SplineFilter_gapfill.R")
library(ggplot2)
#remove the filters--Take DE-Lnf(e.g.EVI) as the example
df.filter<-c()
site.names<-unique(df.VIs$sitename)
for(i in 1:length(site.names)){
  df.proc<-df.VIs %>% 
    filter(sitename==site.names[i])
  #start from non-NA;
  pos.minNA<-min(which(!is.na(df.proc$EVI)))
  df.proc<-df.proc[pos.minNA:nrow(df.proc),]
  df.proc.VIs<-df.proc[,c("EVI","NDVI","NIRv","kNDVI")]
  pos_rm<-match(c("EVI","NDVI","NIRv","kNDVI"),names(df.proc))
  VIs.rm_outliers<-as.data.frame(apply(df.proc.VIs,2,SplineFilter_thengap_f))
  #
  df.proc.new<-cbind(df.proc[,-pos_rm],VIs.rm_outliers)
  #compare the ori ts and filtered ts:
  # gg_EVI<-ggplot()+
  #   geom_point(data=df.proc,aes(x=Date,y=EVI,col="ori"))+
  #   geom_point(data=df.proc.new,aes(x=Date,y=EVI,col="new"))
  # gg_NDVI<-ggplot()+
  #   geom_point(data=df.proc,aes(x=Date,y=NDVI,col="ori"))+
  #   geom_point(data=df.proc.new,aes(x=Date,y=NDVI,col="new"))
  # gg_NIRv<-ggplot()+
  #     geom_point(data=df.proc,aes(x=Date,y=NIRv,col="ori"))+
  #     geom_point(data=df.proc.new,aes(x=Date,y=NIRv,col="new"))
  # gg_kNDVI<-ggplot()+
  #     geom_point(data=df.proc,aes(x=Date,y=kNDVI,col="ori"))+
  #     geom_point(data=df.proc.new,aes(x=Date,y=kNDVI,col="new"))
  df.filter<-rbind(df.filter,df.proc.new)
}

#-----2c).select the site-year data that eligible for phenology extraction:-----
#--select the years that at least have half year data
df.filter_sel<-c()
for(i in 1:length(site.names)){
  df.proc<-df.filter %>% 
    mutate(Year=year(Date),
           Month=month(Date))%>%
    filter(sitename==site.names[i])
  avai.Years<-unique(df.proc$Year)
  #
  df.proc.new<-c()
  for(j in 1:length(avai.Years)){
    temp<-df.proc %>% 
      filter(Year==avai.Years[j])
    avai.data.N<-length(!is.na(temp$EVI))
    #if half year data is not available then do not use this year data
    if(c(avai.data.N/365)<0.5) break
    #test:
    ggplot(data=temp,aes(x=Date,y=EVI))+
      geom_point()
    df.proc.new<-rbind(df.proc.new,temp)
  }
  df.filter_sel<-rbind(df.filter_sel,df.proc.new)
}

##-----2d) add PFTs information--------------
site_info<-read.csv(file="./data-raw/raw_data/SiteInfo.csv")
names(site_info)<-c("sitename","nYear","nYearMODIS","PFT","Latitude","Longitude" )
site_info<-site_info %>%
  dplyr::select(sitename,PFT,Latitude,Longitude)
#
df.filter_sel<-left_join(df.filter_sel,site_info,by="sitename")


#-------2e) start to extract the phenology
source("./R/pheno_extraction_fun2.R")
Phenos_all<-c()
for(i in 1:length(site.names)){
  df.proc<-df.filter_sel %>%
    filter(sitename==site.names[i]) %>%
    mutate(DoY=lubridate::yday(Date))
  #start to extract the phenology:
  Years<-unique(df.proc$Year)
  #
  Phenos_merge<-c()
  site.name<-site.names[i]
  for(j in 1:length(Years)){
    Pheno_results<-SplinePheno_extraction(df.proc,site.name,"EVI",do_norm = FALSE,Years[j])
    #time-series and phenological dates
    SplinePheno_metrics_plot(Pheno_results,site.name,"EVI",Years[j])
    #tidy the phenos for each year
    phenos_temp<-Pheno_results$Pheno_sum$pheno
    phenos_temp$Year<-Years[j]
    Phenos_merge<-rbind(Phenos_merge,phenos_temp)
  }
  Phenos_merge$sitename<-rep(site.name,nrow(Phenos_merge))
  Phenos_all<-bind_rows(Phenos_all,Phenos_merge)
}

#------------------------------------
#3) save the data
#-----------------------------------
save(Phenos_all,file = paste0("./data/data_prep/Extracted_Phenos/df_VIs_Phenos_for_Lu2022sites.rda"))
