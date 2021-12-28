######################################################
##Aim:To tidy the MODIS VIs from Walther et al., 2021:
#https://bg.copernicus.org/preprints/bg-2021-314/
######################################################
library(dplyr)
#----------------
#(1) sites needs to be check for the VIs
#----------------
load(paste0("./data/flux_meteo_VIs_fromJiangong.rda"))
unique(df$sitename)#the used sites are mainly located in Europe, North America

#--------------
#(2) download the MODIS data from Walther et al., 2021
#-------------
#From MODIS: https://meta.icos-cp.eu/collections/tEAkpU6UduMMONrFyym5-tUW

#--------------
#(3) extract the VIs from download nc files-->using the cutout average 
#-------------
#1).manually selected the nc files from those sites:
#a. sites belong to MODIS_Europe_AT_BE_CH_CZ_DK_HU_IE_IS.zip
#b. sites belong to MODIS_Europe_DE.zip
#c. sites belong to MODIS_Europe_ES_FI_FR.zip
#d. sites belong to MODIS_Europe_IT.zip
#e. sites belong to MODIS_Europe_NL_PL_PT_SE_SJ_SK_UK.zip
#f. sites belong to MODIS_Europe_RU.zip
#g. sites belong to MODIS_Americas_CA.zip
#h. sites belong to MODIS_Americas_AR_BR_GF_GL_PA.zip
#i. sites belong to MODIS_Americas_US_1.zip
#j. sites belong to MODIS_Americas_US_2.zip
#k. sites belong to MODIS_Americas_US_3.zip

#2) open the nc files 
library(ncdf4)
library(raster)
library(lubridate)
nc.files.path<-"D:/data/FallPheno_project/MODIS_VIs_from_Walther2021/nc_files/"
sel_sites<-unique(df$sitename)
df.VIs<-c()
for(i in 1:length(sel_sites)){
  fname<-list.files(path = nc.files.path, pattern=sel_sites[i], full.names = TRUE)
  sitename<-sel_sites[i]
  if (length(fname) > 0){
    print(paste('---> Analyzing ', sitename))
    nc <- nc_open(fname)
    
    ##a. extract the variables from the nc file
    #Date-->since 1970-01-01
    N<-nc$dim$time$len
    Date<-as.Date(nc$dim$time$vals,origin="1970-01-01")
    
    #VIs-->BLUE,GREEN,RED,NIR,EVI,NDVI,NIRv,generalized NDVI(kNDVI)
    b_blue<-ncvar_get(nc,"BLUE")
    b_green<-ncvar_get(nc,"GREEN")
    b_red<-ncvar_get(nc,"RED")
    b_nir<-ncvar_get(nc,"NIR")
    EVI<-ncvar_get(nc,"EVI")
    NDVI<-ncvar_get(nc,"NDVI")
    NIRv<-ncvar_get(nc,"NIRv")
    kNDVI<-ncvar_get(nc,"kNDVI")
    
    ##b.create the df
    #produce the new datasets
    df.temp<-data.frame(sitename=rep(sitename,N),Date=Date,b_blue=b_blue,b_green=b_green,b_red=b_red,
                         b_nir=b_nir,EVI=EVI,NDVI=NDVI,NIRv=NIRv,kNDVI=kNDVI)
    df.VIs<-rbind(df.VIs,df.temp)
  }
}
##testing
library(ggplot2)
df.VIs %>% filter(sitename=="DE-Hai") %>%
  ggplot(.,aes(x=Date,y=NDVI))+
  geom_point()
#comparing different data sources
ggplot()+
  geom_point(data = df[df$sitename=="DE-Hai",],aes(x=TIMESTAMP,y=NDVI,col="Jiangong"))+
  geom_point(data = df.VIs[df.VIs$sitename=="DE-Hai",],aes(x=Date,y=NDVI,col="Sophia"))+
  annotate(geom = "text",as.Date("2018-01-01"),y=0.25,label="DE-Hai")
  
#save the MODIS VIs
# save(df.VIs,file = paste0("./data/MODIS_VIs_fromWalther2021.rda"))

#--------------
#(4)Merge the VIs(using data from Walther et al., 2021), and GPP and Meteo(from Jiangong)
#-------------
df<-df %>% mutate(Date=TIMESTAMP) %>%
  mutate(TIMESTAMP=NULL,NDVI=NULL,EVI=NULL,EVI_2=NULL,NIRv=NULL)
df.VIs<-df.VIs %>%
  mutate(Year=year(Date))
df.final <- left_join(df, df.VIs,by=c("sitename","Year","Date"))
#save the final data
save(df.final,file = paste0("./data/df_GPP_Meteo_andVIs.rda"))
