############################################################
#Aim: re-tidy the variables(like in script 1_Data...) 
#to explore the relationship between EOS and its drivers
############################################################
library(dplyr)
library(tidyverse)
#--------------------
#1.load the data
#--------------------
load(paste0("./data/data_prep/Extracted_Phenos/df_VIs_Phenos_updated.rda")) #phenology data
load(paste0("./data/data_prep/df_GPP_Meteo_andVIs.rda")) #datasets(GPP,Meteo,and VIs)

#-------------------------
#2.merge the data and select the site-years that available for the analysis
#------------------------
df.merge<-left_join(df.final,Phenos_final,by=c("sitename","Year","PFT"))
#only select the site-years that reasonable for the analysis:
#first of all, remove the site-year that sos and eos are not reliable:
df.sel<-df.merge %>%
  filter(!is.na(trs_sos25)&!is.na(trs_eos25))%>% #sos25 and eos25 are both available#sos25 and eos25 are both available
  filter(trs_sos25>=60 & trs_sos25<=173)%>% #sos should between March and summer solstice
  filter(trs_eos90>173) %>% #eos should later than summer solstice
  mutate(summer_sol=173,
         GSL_sos25sol=summer_sol - trs_sos25,GSL_solpeak=pop - summer_sol,
         GSL_peakeos90=trs_eos90 - pop,GSL_eos90eos50=trs_eos50 - trs_eos90,
         GSL_eos50eos25=trs_eos25 - trs_eos50,  ##for Period A (A1-A5)
         GSL_sos25peak=pop - trs_sos25,GSL_sos25eos90=trs_eos90 - trs_sos25,
         GSL_sos25eos50=trs_eos50 - trs_sos25,
         GSL_sos25eos25=trs_eos25 - trs_sos25  ##for Period B(B1-B4)
         )

#-------------------------------------------
#3.calculate the variables
#mean of (Meterological variables):
#temperature (Ta,Tday,Tnight,Tmax,Tmin)
#SW_IN,VPD
#mean and cum GPP_NT,GPP_DT,NEP
#in different periods in Period A and Period B
#-------------------------------------------
df.proc<-df.sel %>%
  select(sitename,Date,Year,Month,PFT,
         TA_F,TA_F_DAY,TA_F_NIGHT,TA_min,TA_max,
         SW_IN_F,VPD_F,
         GPP_NT_VUT_REF,GPP_DT_VUT_REF,NEE_VUT_USTAR50,
         starts_with("trs_"),summer_sol,pop,
         starts_with("GSL_")
         ) %>%
  mutate(DoY=lubridate::yday(Date))

###
agg_fun<-function(df,sel_period){
  # df<-df.proc
  # sel_period<-"sos25sol"
  
  #specify the period:
  ##Period A:
  if(sel_period=="sos25sol"){
    df.run<-df %>%
      filter(DoY>=trs_sos25 & DoY<=summer_sol)
  }
  if(sel_period=="solpeak"){
    df.run<-df %>%
      filter(DoY>=summer_sol & DoY<=pop)
  }
  if(sel_period=="peakeos90"){
    df.run<-df %>%
      filter(DoY>=pop & DoY<=trs_eos90)
  }
  if(sel_period=="eos90eos50"){
    df.run<-df %>%
      filter(DoY>=trs_eos90 & DoY<=trs_eos50)
  }
  if(sel_period=="eos50eos25"){
    df.run<-df %>%
      filter(DoY>=trs_eos50 & DoY<=trs_eos25)
  }
  ##Period B:
  if(sel_period=="sos25peak"){
    df.run<-df %>%
      filter(DoY>=trs_sos25 & DoY<=pop)
  }
  if(sel_period=="sos25eos90"){
    df.run<-df %>%
      filter(DoY>=trs_sos25 & DoY<=trs_eos90)
  }
  if(sel_period=="sos25eos50"){
    df.run<-df %>%
      filter(DoY>=trs_sos25 & DoY<=trs_eos50)
  }
  if(sel_period=="sos25eos25"){
    df.run<-df %>%
      filter(DoY>=trs_sos25 & DoY<=trs_eos25)
  }
  
  #calculate the stats:
  df.stats<-df.run %>%
    group_by(sitename,Year) %>%
    dplyr::summarise(Tmean=mean(TA_F,na.rm=T),
        Tmin=mean(TA_min,na.rm=T),Tmax=mean(TA_max,na.rm=T),
        Tday=mean(TA_F_DAY,na.rm=T),Tnight=mean(TA_F_NIGHT),
        SW_IN=mean(SW_IN_F,na.rm=T),VPD=mean(VPD_F,na.rm=T),
        GPP_NT_mean=mean(GPP_NT_VUT_REF,na.rm=T),
        GPP_DT_mean=mean(GPP_DT_VUT_REF,na.rm=T),
        NPP_mean=c(-1)*mean(NEE_VUT_USTAR50,na.rm=T)
        )
  df.phenos<-df.run %>%
    select(sitename,Year,trs_sos25:GSL_sos25eos25)%>%
    group_by(sitename,Year) %>%
    dplyr::summarise_at(vars(trs_sos25:GSL_sos25eos25),mean,na.rm=T)
  #merge the datasets and calculate the GPP_cum
  df.out<-left_join(df.stats,df.phenos,by=c("sitename","Year"))
  # pos_GSL=grep(sel_period,names(df.out))
  GSL_name<-paste0("GSL_",sel_period)
  vars<-c("GPP_NT_mean","GPP_DT_mean","NPP_mean",GSL_name)
  df.out<-df.out %>%
    mutate(GPP_NT_cum=.data[[vars[[1]]]]*.data[[vars[[4]]]],
           GPP_DT_cum=.data[[vars[[2]]]]*.data[[vars[[4]]]],
           NEP_cum=.data[[vars[[3]]]]*.data[[vars[[4]]]]
           )
  #
  return(df.out)
}

##summary for different periods:
df.sos25sol<-agg_fun(df.proc,"sos25sol")


df.solpeak<-agg_fun(df.proc,"solpeak")
