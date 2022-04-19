##############################
#Aim: to tidy and merge the data the shared data from Lu et al., 2022:
#https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16104
##############################
library(dplyr)
library()
#------------
proc.path<-"./data-raw/raw_data/data_from_Luetal2022/FLX_ALL_phe/"
setwd(proc.path)
files<-list.files()
#site names:
sites<-substr(files,1,6)
#open the csv files and merge the data:
df_merge<-c()
for (i in 1:length(sites)) {
  temp<-read.csv(files[i])
  temp_add<-data.frame(sitename=rep(sites[i],nrow(temp)),temp)
  df_merge<-bind_rows(df_merge,temp_add)  
  rm(temp)
}
#save the data
save.path<-"./data/data_prep/Add_analysis_based_Lu2022/"
setwd("D:/Github/FallPheno")
save(df_merge,file = paste0(save.path,"Data_from_Lu2022.RDA"))
