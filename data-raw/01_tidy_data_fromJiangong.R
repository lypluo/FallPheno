#################################################################
#tidy the Fluxes, Meterological variables and VIs from Jiangong
#################################################################
Sys.setenv(tz="UTC")
library(dplyr)
#---------
#1) load the site info
#---------
site_info<-read.csv(file="./data-raw/raw_data/SiteInfo.csv")
par(mfrow=c(1,2))
hist(site_info$nYearMODIS)
hist(site_info$nYear)
#select the sites have more year 6,8,10 years fluxes and modis data
# sel_sites1<-site_info %>% 
#   filter(nYear > 6 & nYearMODIS >6) #88
# sel_sites2<-site_info %>% 
#   filter(nYear > 8 & nYearMODIS >8) #68
# sel_sites3<-site_info %>% 
#   filter(nYear > 10 & nYearMODIS >10) #55

### Select the nYear and nYearMODIS >10
sel_sites<-site_info %>%
  filter(nYear > 10 & nYearMODIS >10) 
sel_sites %>%
  group_by(PFT) %>%
  summarise(nsites=length(PFT))
### First focus on the Main natural PFTs-->DBF, ENF,GRA,MF
final_sites<-sel_sites %>%
  filter(PFT %in% c("DBF","ENF","GRA","MF")) #37sites
#save the sel site info
write.csv(final_sites,file = paste0("./data/data_prep/SiteInfo_sel.csv"))

#--------------------------
#(2)load the flux,VIs and meterolgoical data
#---------------------------
#path for the datasets
df.path<-"D:/data/FallPheno_project/orgin_data_from_Jiangong/Data_Liu/"
#
df<-c()
for(i in 1:nrow(final_sites)){
  site.name<-final_sites$Site[i]
  df.temp<-read.csv(paste0(df.path,site.name,"_DD_FLUXNET_VI.csv"))
  df.temp<- df.temp %>% dplyr::select(-c(YYDD,Day)) #remove the variable Day
  #
  N<-nrow(df.temp)
  df.new<- data.frame(sitename=rep(site.name,N),Flux_nYear=rep(final_sites$nYear[i],N),
                      MODIS_nYear=rep(final_sites$nYearMODIS[i],N),PFT=rep(final_sites$PFT[i],N))
  df.new<-cbind(df.new,df.temp)
  df<-bind_rows(df,df.new)
}

#convert the formate of the date
df$TIMESTAMP<-as.Date(df$TIMESTAMP)

#-----------------------
#merge site info and data
#-----------------------
final_sites<-final_sites %>%mutate(sitename=Site,Site=NULL)
df<-left_join(df,final_sites)
  
#save the data
save(df,file = paste0("./data/data_prep/EC_data/flux_meteo_VIs_fromJiangong.rda"))





