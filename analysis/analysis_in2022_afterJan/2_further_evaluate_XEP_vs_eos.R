#####################################
#Aim: Further evaluate what's the relationship between GPP/NEP mean or cum vs eos
#####################################
library(dplyr)
#--------------------
#1.load the data ca
#--------------------
load(paste0("./data/df_VIs_Phenos_updated.rda")) #phenology data
load(paste0("./data/df_GPP_Meteo_andVIs.rda")) #datasets(GPP,Meteo,and VIs)

#--------------------
#2.select the data could be used for analysis -->for GPP, Meteo, and VIs dataset
#--------------------
#select the data for analysis
df.analysis<-c()
for(i in 1:nrow(Phenos_final)){
  Pheno.temp<-Phenos_final[i,]
  df.temp<-df.final %>%
    filter(sitename==Pheno.temp$sitename & Year==Pheno.temp$Year)
  df.analysis<-rbind(df.analysis,df.temp)
}
#b.check the univariate correlation between phenophases:
library(corrplot)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#--------------------
#3.calculate the GPP mean, NEE daily mean, Cum GPP/NEE during the growing season in each site-year
#here calculating those variables in different time period:
#Apart from:[sos25,solstice],[sos25,eos90],[sos25,eos50],[sos25,eos25]
#we also add the analysis for following period on Feb,2022:[solstice,eos90],
#[eos90,eos50],[eos50,eos25]
#--------------------
library(lubridate)
df.sum<-c()
#first of all, remove the site-year that sos and eos are not reliable:
Phenos_final_sel<-Phenos_final%>%
  filter(!is.na(trs_sos25)&!is.na(trs_eos25))%>% #sos25 and eos25 are both available
  filter(trs_sos25>=60 & trs_sos25<=173)%>% #sos should between March and summer solstice
  filter(trs_eos90>173) %>%#eos should later than summer solstice
  mutate(summer_sol=173,
         GSL_sos25sol=summer_sol - trs_sos25,GSL_sos25eos90=trs_eos90 - trs_sos25,
         GSL_sos25eos50=trs_eos50 - trs_sos25,GSL=trs_eos25 - trs_sos25,
         GSL_soleos90=trs_eos90 - summer_sol,GSL_eos90eos50= trs_eos50 - trs_eos90,#add more period length
         GSL_eos50eos25=trs_eos25 - trs_eos50
         )

for(i in 1:nrow(Phenos_final_sel)){
  Pheno.temp<-Phenos_final_sel[i,]
  #
  Phenos.sel<-Pheno.temp %>%
    select(sitename,Year,PFT,starts_with("trs_sos"),pop,starts_with("trs_eos"),
           summer_sol,starts_with("GSL"),
           maxline,baseline,ampl_sos,ampl_eos,prr,psr,Gslope,Dslope)
  df.temp<-df.analysis %>%
    filter(sitename==Pheno.temp$sitename & Year==Pheno.temp$Year) %>%
    mutate(DoY=lubridate::yday(Date)) %>%
    filter(DoY>=Pheno.temp$trs_sos25 & DoY<=Pheno.temp$trs_eos25)
  #merge:
  df.temp<-left_join(df.temp,Phenos.sel,by=c("sitename","Year"))
  ##here re
  df.sum.phenos<-df.temp %>%
    filter(DoY>=trs_sos25 & DoY<=trs_eos25)%>%   #focus on the growing season: [sos25,eos25] 
    summarise(sitename=unique(sitename),Year=unique(Year),
              sos25=mean(trs_sos25,na.rm=T),summer_solstice=173,pop=mean(pop,na.rm=T), #set summer solstice equals 6-22
              eos90=mean(trs_eos90,na.rm=T),eos50=mean(trs_eos50,na.rm=T),eos25=mean(trs_eos25,na.rm=T),
              GSL_sos25sol=mean(GSL_sos25sol,na.rm=T),GSL_sos25eos90=mean(GSL_sos25eos90,na.rm=T),
              GSL_sos25eos50=mean(GSL_sos25eos50,na.rm=T),GSL=mean(GSL,na.rm=T),
              GSL_soleos90=mean(GSL_soleos90,na.rm=T),GSL_eos90eos50=mean(GSL_eos90eos50,na.rm=T), #add more period length
              GSL_eos50eos25=mean(GSL_eos50eos25,na.rm=T)
              )
  ##------for different periods-->sum GPP/NEP mean
  df.sum.sos25sol<-df.temp %>%
    filter(DoY>=trs_sos25 & DoY<=173) %>%
    summarise(GPP_NT_GSL25sol_mean=mean(GPP_NT_VUT_REF,na.rm=T), #-->GPP mean for different period
              GPP_NT_GSL25sol_mean=mean(GPP_NT_VUT_REF,na.rm=T),   #unit:g C m-2 d-1
              GPP_DT_GSL25sol_mean=mean(GPP_DT_VUT_REF,na.rm=T),
              NEP_GSL25sol_mean=c(-mean(NEE_VUT_USTAR50,na.rm=T)))
  
  df.sum.sos25eos90<-df.temp %>%
    filter(DoY>=trs_sos25 & DoY<=trs_eos90) %>%
    summarise(GPP_NT_GSL2590_mean=mean(GPP_NT_VUT_REF,na.rm=T), #-->GPP mean for different period
              GPP_NT_GSL2590_mean=mean(GPP_NT_VUT_REF,na.rm=T),   #unit:g C m-2 d-1
              GPP_DT_GSL2590_mean=mean(GPP_DT_VUT_REF,na.rm=T),
              NEP_GSL2590_mean=c(-mean(NEE_VUT_USTAR50,na.rm=T)))
  
  df.sum.sos25eos50<-df.temp %>%
    filter(DoY>=trs_sos25 & DoY<=trs_eos50) %>%
    summarise(GPP_NT_GSL2550_mean=mean(GPP_NT_VUT_REF,na.rm=T), #-->GPP mean for different period
              GPP_NT_GSL2550_mean=mean(GPP_NT_VUT_REF,na.rm=T),   #unit:g C m-2 d-1
              GPP_DT_GSL2550_mean=mean(GPP_DT_VUT_REF,na.rm=T),
              NEP_GSL2550_mean=c(-mean(NEE_VUT_USTAR50,na.rm=T)))
  df.sum.sos25eos25<-df.temp %>%
    filter(DoY>=trs_sos25 & DoY<=trs_eos25) %>%
    summarise(GPP_NT_GSL_mean=mean(GPP_NT_VUT_REF,na.rm=T), #-->GPP mean for different period
              GPP_NT_GSL_mean=mean(GPP_NT_VUT_REF,na.rm=T),   #unit:g C m-2 d-1
              GPP_DT_GSL_mean=mean(GPP_DT_VUT_REF,na.rm=T),
              NEP_GSL_mean=c(-mean(NEE_VUT_USTAR50,na.rm=T)))
  ##adding more periods data:
  df.sum.soleos90<-df.temp %>%
    filter(DoY>=173 & DoY<=trs_eos90) %>%
    summarise(GPP_NT_GSLsol90_mean=mean(GPP_NT_VUT_REF,na.rm=T), #-->GPP mean for different period
              GPP_NT_GSLsol90_mean=mean(GPP_NT_VUT_REF,na.rm=T),   #unit:g C m-2 d-1
              GPP_DT_GSLsol90_mean=mean(GPP_DT_VUT_REF,na.rm=T),
              NEP_GSLsol90_mean=c(-mean(NEE_VUT_USTAR50,na.rm=T)))
  df.sum.eos9050<-df.temp %>%
    filter(DoY>=trs_eos90 & DoY<=trs_eos50) %>%
    summarise(GPP_NT_GSL9050_mean=mean(GPP_NT_VUT_REF,na.rm=T), #-->GPP mean for different period
              GPP_NT_GSL9050_mean=mean(GPP_NT_VUT_REF,na.rm=T),   #unit:g C m-2 d-1
              GPP_DT_GSL9050_mean=mean(GPP_DT_VUT_REF,na.rm=T),
              NEP_GSL9050_mean=c(-mean(NEE_VUT_USTAR50,na.rm=T)))
  df.sum.eos5025<-df.temp %>%
    filter(DoY>=trs_eos50 & DoY<=trs_eos25) %>%
    summarise(GPP_NT_GSL5025_mean=mean(GPP_NT_VUT_REF,na.rm=T), #-->GPP mean for different period
              GPP_NT_GSL5025_mean=mean(GPP_NT_VUT_REF,na.rm=T),   #unit:g C m-2 d-1
              GPP_DT_GSL5025_mean=mean(GPP_DT_VUT_REF,na.rm=T),
              NEP_GSL5025_mean=c(-mean(NEE_VUT_USTAR50,na.rm=T)))
  
  ##---merge together:
  df.sum.temp<-cbind(df.sum.phenos,df.sum.sos25sol,df.sum.sos25eos90,
                     df.sum.sos25eos50,df.sum.sos25eos25,
                     df.sum.soleos90,df.sum.eos9050,df.sum.eos5025 ##add more period
                     )
  df.tidy<-df.sum.temp %>%      #calculate the cum data
    mutate(GPP_NT_GSL25sol_cum = GPP_NT_GSL25sol_mean*GSL_sos25sol,
           GPP_DT_GSL25sol_cum=GPP_DT_GSL25sol_mean*GSL_sos25sol,
           NEP_GSL25sol_cum=NEP_GSL_mean*GSL_sos25sol, #
           GPP_NT_GSL2590_cum = GPP_NT_GSL2590_mean*GSL_sos25eos90,
           GPP_DT_GSL2590_cum=GPP_DT_GSL2590_mean*GSL_sos25eos90,
           NEP_GSL2590_cum=NEP_GSL_mean *GSL_sos25eos90,#
           GPP_NT_GSL2550_cum = GPP_NT_GSL2550_mean*GSL_sos25eos50,
           GPP_DT_GSL2550_cum=GPP_DT_GSL2550_mean*GSL_sos25eos50,
           NEP_GSL2550_cum=NEP_GSL_mean*GSL_sos25eos50,#
           GPP_NT_GSL_cum = GPP_NT_GSL_mean*GSL,
           GPP_DT_GSL_cum=GPP_DT_GSL_mean*GSL,
           NEP_GSL_cum=NEP_GSL_mean*GSL,
           GPP_NT_GSLsol90_cum = GPP_NT_GSLsol90_mean*GSL_soleos90, #add more period data
           GPP_DT_GSLsol90_cum = GPP_DT_GSLsol90_mean*GSL_soleos90,
           NEP_GSLsol90_cum = NEP_GSLsol90_mean*GSL_soleos90,
           GPP_NT_GSL9050_cum = GPP_NT_GSL9050_mean*GSL_eos90eos50, #add more period data
           GPP_DT_GSL9050_cum = GPP_DT_GSL9050_mean*GSL_eos90eos50,
           NEP_GSL9050_cum = NEP_GSL9050_mean*GSL_eos90eos50,
           GPP_NT_GSL5025_cum = GPP_NT_GSL5025_mean*GSL_eos50eos25, #add more period data
           GPP_DT_GSL5025_cum = GPP_DT_GSL5025_mean*GSL_eos50eos25,
           NEP_GSL5025_cum = NEP_GSL5025_mean*GSL_eos50eos25
    )
  df.sum<-rbind(df.sum,df.tidy)
}

#----
#!!filter the outliers-->remove the site-years has negative GPP_NT_mean or GPP_DT_mean
#----
df.outlier<-df.sum %>%
  filter(GPP_NT_GSL_mean<=0 |GPP_DT_GSL_mean<=0)
#outliers-->US-WCr: 2010; US-Syv:2008;GL-Zah:2006
df.merge.new<-df.sum %>%
  filter(GPP_NT_GSL_mean>0 & GPP_DT_GSL_mean>0)  ## also remove when GPP value == NA

#site:GL-ZaH,US-Wkg.
#----
#remove the sites that have very low GPP mean
#----
df.merge.final<-df.merge.new %>%
  filter(GPP_NT_GSL_mean > 3 & GPP_DT_GSL_mean >3)
plot(x=df.merge.final$GPP_NT_GSL_mean,y=df.merge.final$eos25)
plot(x=df.merge.final$GPP_DT_GSL_mean,y=df.merge.final$eos25)
plot(x=df.merge.final$GPP_NT_GSL_cum,y=df.merge.final$eos25)
plot(x=df.merge.final$GPP_DT_GSL_cum,y=df.merge.final$eos25)
#adding more
plot(x=df.merge.final$GPP_NT_GSLsol90_mean,y=df.merge.final$eos90)
plot(x=df.merge.final$GPP_NT_GSL9050_mean,y=df.merge.final$eos50)
plot(x=df.merge.final$GPP_NT_GSL5025_mean,y=df.merge.final$eos25)
plot(x=df.merge.final$GPP_NT_GSLsol90_cum,y=df.merge.final$eos90)
plot(x=df.merge.final$GPP_NT_GSL9050_cum,y=df.merge.final$eos50)
plot(x=df.merge.final$GPP_NT_GSL5025_cum,y=df.merge.final$eos25)

#----------------------------
#(5) evaluate the realtionship between gpp/nep vs eos for different periods
#----------------------------
library(ggplot2)
library(ggpubr)
#------------------------
##I.using all the original values:
#-----------------------
df.new<-df.merge.final
##1.correlation plot-->sos,eos
#1a. for phenophase, daily average GPP/NEP, and cum GPP/NEP
# df.new_sel1<-df.new[,c("sos25","eos90","eos50","eos25","pop",
#                        "GPP_NT_GSL_mean","GPP_DT_GSL_mean","NEP_GSL_mean",
#                        "GPP_NT_GSL25sol_mean","GPP_DT_GSL25sol_mean","NEP_GSL25sol_mean",
#                        "GPP_NT_GSL2550_mean","GPP_DT_GSL2550_mean","NEP_GSL2550_mean",
#                        "GPP_NT_GSL2590_mean","GPP_DT_GSL2590_mean","NEP_GSL2590_mean",
#                        "GPP_NT_GSL_cum","GPP_DT_GSL_cum","NEP_GSL_cum",
#                        "GPP_NT_GSL25sol_cum","GPP_DT_GSL25sol_cum","NEP_GSL25sol_cum",
#                        "GPP_NT_GSL2550_cum","GPP_DT_GSL2550_cum","NEP_GSL2550_cum",
#                        "GPP_NT_GSL2590_cum","GPP_DT_GSL2590_cum","NEP_GSL2590_cum")]
df.new_sel1<-df.new[,c("sos25","eos90","eos50","eos25","pop",
    "GPP_NT_GSL25sol_mean","GPP_DT_GSL25sol_mean","NEP_GSL25sol_mean",
    "GPP_NT_GSL25sol_cum","GPP_DT_GSL25sol_cum","NEP_GSL25sol_cum")]
df.new_sel2<-df.new[,c("sos25","eos90","eos50","eos25","pop",
                       "GPP_NT_GSL2590_mean","GPP_DT_GSL2590_mean","NEP_GSL2590_mean",
                       "GPP_NT_GSL2590_cum","GPP_DT_GSL2590_cum","NEP_GSL2590_cum")]
df.new_sel3<-df.new[,c("sos25","eos90","eos50","eos25","pop",
    "GPP_NT_GSL2550_mean","GPP_DT_GSL2550_mean","NEP_GSL2550_mean",
    "GPP_NT_GSL2550_cum","GPP_DT_GSL2550_cum","NEP_GSL2550_cum")]
df.new_sel4<-df.new[,c("sos25","eos90","eos50","eos25","pop",
                       "GPP_NT_GSL_mean","GPP_DT_GSL_mean","NEP_GSL_mean",
                       "GPP_NT_GSL_cum","GPP_DT_GSL_cum","NEP_GSL_cum")]
##adding more the periods:
df.new_sel5<-df.new[,c("sos25","eos90","eos50","eos25","pop",
                       "GPP_NT_GSLsol90_mean","GPP_DT_GSLsol90_mean","NEP_GSLsol90_mean",
                       "GPP_NT_GSLsol90_cum","GPP_DT_GSLsol90_cum","NEP_GSLsol90_cum")]
df.new_sel6<-df.new[,c("sos25","eos90","eos50","eos25","pop",
                       "GPP_NT_GSL9050_mean","GPP_DT_GSL9050_mean","NEP_GSL9050_mean",
                       "GPP_NT_GSL9050_cum","GPP_DT_GSL9050_cum","NEP_GSL9050_cum")]
df.new_sel7<-df.new[,c("sos25","eos90","eos50","eos25","pop",
                       "GPP_NT_GSL5025_mean","GPP_DT_GSL5025_mean","NEP_GSL5025_mean",
                       "GPP_NT_GSL5025_cum","GPP_DT_GSL5025_cum","NEP_GSL5025_cum")]

#------------correlation matrixes-------
M1<-cor(df.new_sel1,use = "complete.obs")
p.mat.M1<-cor.mtest(df.new_sel1)
M2<-cor(df.new_sel2,use = "complete.obs")
p.mat.M2<-cor.mtest(df.new_sel2)
M3<-cor(df.new_sel3,use = "complete.obs")
p.mat.M3<-cor.mtest(df.new_sel3)
M4<-cor(df.new_sel4,use = "complete.obs")
p.mat.M4<-cor.mtest(df.new_sel4)
M5<-cor(df.new_sel5,use = "complete.obs") #adding more
p.mat.M5<-cor.mtest(df.new_sel5)
M6<-cor(df.new_sel6,use = "complete.obs")
p.mat.M6<-cor.mtest(df.new_sel6)
M7<-cor(df.new_sel7,use = "complete.obs")
p.mat.M7<-cor.mtest(df.new_sel7)

#refer the code here: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogramvo
png(filename = paste("./fig/Results_updated/using_filtered_data/uni_corrlation_PhenovsGPP_additional.png"),
    width = 1024,height = 1024)
# layout(matrix(1,2,3,4),2,2)
par(mfrow=c(4,2),mai=c(2,2,4,2))
corrplot(M1, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M1,sig.level = 0.01,addCoef.col = "black",insig = "blank")
# par(fig=c(0.5,1,0.5,1),new=T)
corrplot(M2, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M2,sig.level = 0.01,addCoef.col = "black",insig = "blank")
# par(fig=c(0,0.5,0,0.5),new=T)
corrplot(M3, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M3,sig.level = 0.01,addCoef.col = "black",insig = "blank")
# par(fig=c(0.5,1,0,0.5),new=T)
corrplot(M4, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M4,sig.level = 0.01,addCoef.col = "black",insig = "blank")
##adding more:
corrplot(M5, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M5,sig.level = 0.01,addCoef.col = "black",insig = "blank")
corrplot(M6, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M6,sig.level = 0.01,addCoef.col = "black",insig = "blank")
corrplot(M7, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M7,sig.level = 0.01,addCoef.col = "black",insig = "blank")
dev.off()



####2.regression plot
#(1) mean and cum GPP/NEP vs eos
##########################
#mean GPP/NEP vs eos
##########################
#A.eos90 ~ GPP/NEP in [sos25,solstice]
p_eos90_vs_GPP_NT_GSL25sol<-ggscatter(df.new,
  x="GPP_NT_GSL25sol_mean",y="eos90",add = "reg.line")+
    stat_smooth(method = "lm",formula = y ~ x)+
    stat_cor(label.x = 1, label.y = 185,col="blue") +
    stat_regline_equation(label.x = 1, label.y =172,col="blue" )
lm_temp<-lm(data=df.new,
            eos90 ~ GPP_NT_GSL25sol_mean)
summary(lm_temp)
p_eos90_vs_GPP_DT_GSL25sol<-ggscatter(df.new,
                           x="GPP_DT_GSL25sol_mean",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 1, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 1, label.y = 172,col="blue")
p_eos90_vs_NEP_GSL25sol<-ggscatter(df.new,
                           x="NEP_GSL25sol_mean",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")

#B.eos90 ~ GPP/NEP in [sos25,eos90]
p_eos90_vs_GPP_NT_GSL2590<-ggscatter(df.new,
                                      x="GPP_NT_GSL2590_mean",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue" )
lm_temp<-lm(data=df.new,
            eos90 ~ GPP_NT_GSL2590_mean)
summary(lm_temp)
p_eos90_vs_GPP_DT_GSL2590<-ggscatter(df.new,
                                      x="GPP_DT_GSL2590_mean",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_eos90_vs_NEP_GSL2590<-ggscatter(df.new,
                                   x="NEP_GSL2590_mean",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")
#C.eos50 ~ GPP/NEP in [sos25,eos50]
p_eos50_vs_GPP_NT_GSL2550<-ggscatter(df.new,
                                     x="GPP_NT_GSL2550_mean",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue" )
lm_temp<-lm(data=df.new,
            eos50 ~ GPP_NT_GSL2550_mean)
summary(lm_temp)
p_eos50_vs_GPP_DT_GSL2550<-ggscatter(df.new,
                                     x="GPP_DT_GSL2550_mean",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_eos50_vs_NEP_GSL2550<-ggscatter(df.new,
                                  x="NEP_GSL2550_mean",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")
#D.eos25 ~ GPP/NEP in [sos25,eos25]
p_eos25_vs_GPP_NT_GSL<-ggscatter(df.new,
                                     x="GPP_NT_GSL_mean",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
lm_temp<-lm(data=df.new,
            eos25 ~ GPP_NT_GSL_mean)
summary(lm_temp)
p_eos25_vs_GPP_DT_GSL<-ggscatter(df.new,
                                     x="GPP_DT_GSL_mean",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_eos25_vs_NEP_GSL<-ggscatter(df.new,
                                  x="NEP_GSL_mean",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")
#-------adding more periods data:
#E.eos25 ~ GPP/NEP in [sosltice,eos90]
p_eos90_vs_GPP_NT_GSLsol90<-ggscatter(df.new,
                                 x="GPP_NT_GSLsol90_mean",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
lm_temp<-lm(data=df.new,
            eos90 ~ GPP_NT_GSLsol90_mean)
summary(lm_temp)
p_eos90_vs_GPP_DT_GSLsol90<-ggscatter(df.new,
                                    x="GPP_DT_GSLsol90_mean",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
p_eos90_vs_NEP_GSLsol90<-ggscatter(df.new,
                              x="NEP_GSLsol90_mean",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")

#F.eos50 ~ GPP/NEP in [eos90,eos50]
p_eos50_vs_GPP_NT_GSL9050<-ggscatter(df.new,
                                    x="GPP_NT_GSL9050_mean",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
p_eos50_vs_GPP_DT_GSL9050<-ggscatter(df.new,
                                    x="GPP_DT_GSL9050_mean",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
p_eos50_vs_NEP_GSL9050<-ggscatter(df.new,
                                 x="NEP_GSL9050_mean",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")

#G.eos25 ~ GPP/NEP in [eos50,eos25]
p_eos25_vs_GPP_NT_GSL5025<-ggscatter(df.new,
                                   x="GPP_NT_GSL5025_mean",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
p_eos25_vs_GPP_DT_GSL5025<-ggscatter(df.new,
                                   x="GPP_DT_GSL5025_mean",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
p_eos25_vs_NEP_GSL5025<-ggscatter(df.new,
                                x="NEP_GSL5025_mean",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")

##########################
#cum GPP/NEP vs eos
##########################
p_cum_eos90_vs_GPP_NT_GSL25sol<-ggscatter(df.new,
                                      x="GPP_NT_GSL25sol_cum",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 1, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 1, label.y =172,col="blue")
lm_temp<-lm(data=df.new,
            eos90 ~ GPP_NT_GSL25sol_cum)
summary(lm_temp)
p_cum_eos90_vs_GPP_DT_GSL25sol<-ggscatter(df.new,
                                      x="GPP_DT_GSL25sol_cum",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 1, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 1, label.y = 172,col="blue")
p_cum_eos90_vs_NEP_GSL25sol<-ggscatter(df.new,
                                   x="NEP_GSL25sol_cum",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")

#B.eos90 ~ GPP/NEP in [sos25,eos90]
p_cum_eos90_vs_GPP_NT_GSL2590<-ggscatter(df.new,
                                     x="GPP_NT_GSL2590_cum",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
lm_temp<-lm(data=df.new,
            eos90 ~ GPP_NT_GSL2590_cum)
summary(lm_temp)
p_cum_eos90_vs_GPP_DT_GSL2590<-ggscatter(df.new,
                                     x="GPP_DT_GSL2590_cum",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_cum_eos90_vs_NEP_GSL2590<-ggscatter(df.new,
                                  x="NEP_GSL2590_cum",y="eos90",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")
#C.eos50 ~ GPP/NEP in [sos25,eos50]
p_cum_eos50_vs_GPP_NT_GSL2550<-ggscatter(df.new,
                                     x="GPP_NT_GSL2550_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
lm_temp<-lm(data=df.new,
            eos50 ~ GPP_NT_GSL2550_cum)
summary(lm_temp)
p_cum_eos50_vs_GPP_DT_GSL2550<-ggscatter(df.new,
                                     x="GPP_DT_GSL2550_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_cum_eos50_vs_NEP_GSL2550<-ggscatter(df.new,
                                  x="NEP_GSL2550_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")
#D.eos25 ~ GPP/NEP in [sos25,eos25]
p_cum_eos25_vs_GPP_NT_GSL<-ggscatter(df.new,
                                 x="GPP_NT_GSL_cum",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
lm_temp<-lm(data=df.new,
            eos25 ~ GPP_NT_GSL_cum)
summary(lm_temp)
p_cum_eos25_vs_GPP_DT_GSL<-ggscatter(df.new,
                                 x="GPP_DT_GSL_cum",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_cum_eos25_vs_NEP_GSL<-ggscatter(df.new,
                              x="NEP_GSL_cum",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")
#-------------adding more data in other periods:
#E.eos90 ~ GPP/NEP in [solstice,eos90]
p_cum_eos90_vs_GPP_NT_GSLsol90<-ggscatter(df.new,
                                         x="GPP_NT_GSLsol90_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
p_cum_eos90_vs_GPP_DT_GSLsol90<-ggscatter(df.new,
                                         x="GPP_DT_GSLsol90_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_cum_eos90_vs_NEP_GSLsol90<-ggscatter(df.new,
                                      x="NEP_GSLsol90_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")

#F.eos50 ~ GPP/NEP in [eos90,eos50]
p_cum_eos50_vs_GPP_NT_GSL9050<-ggscatter(df.new,
                                             x="GPP_NT_GSL9050_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
p_cum_eos50_vs_GPP_DT_GSL9050<-ggscatter(df.new,
                                             x="GPP_DT_GSL9050_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_cum_eos50_vs_NEP_GSL9050<-ggscatter(df.new,
                                          x="NEP_GSL9050_cum",y="eos50",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")

#G.eos25 ~ GPP/NEP in [eos50,eos25]
p_cum_eos25_vs_GPP_NT_GSL5025<-ggscatter(df.new,
                                         x="GPP_NT_GSL5025_cum",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y =172,col="blue")
p_cum_eos25_vs_GPP_DT_GSL5025<-ggscatter(df.new,
                                         x="GPP_DT_GSL5025_cum",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 5, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 5, label.y = 172,col="blue")
p_cum_eos25_vs_NEP_GSL5025<-ggscatter(df.new,
                                      x="NEP_GSL5025_cum",y="eos25",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 0, label.y = 185,col="blue") +
  stat_regline_equation(label.x = 0, label.y = 172,col="blue")

#merge the plots
library(cowplot)
p_eos90_vs_XEP_GSL25sol<-plot_grid(p_eos90_vs_GPP_NT_GSL25sol,
        p_eos90_vs_GPP_DT_GSL25sol,p_eos90_vs_NEP_GSL25sol,
        p_cum_eos90_vs_GPP_NT_GSL25sol,
        p_cum_eos90_vs_GPP_DT_GSL25sol,p_cum_eos90_vs_NEP_GSL25sol)
p_eos90_vs_XEP_GSL2590<-plot_grid(p_eos90_vs_GPP_NT_GSL2590,
        p_eos90_vs_GPP_DT_GSL2590,p_eos90_vs_NEP_GSL2590,
        p_cum_eos90_vs_GPP_NT_GSL2590,
        p_cum_eos90_vs_GPP_DT_GSL2590,p_cum_eos90_vs_NEP_GSL2590)
p_eos50_vs_XEP_GSL2550<-plot_grid(p_eos50_vs_GPP_NT_GSL2550,
        p_eos50_vs_GPP_DT_GSL2550,p_eos50_vs_NEP_GSL2550,
        p_cum_eos50_vs_GPP_NT_GSL2550,
        p_cum_eos50_vs_GPP_DT_GSL2550,p_cum_eos50_vs_NEP_GSL2550)
p_eos25_vs_XEP_GSL<-plot_grid(p_eos25_vs_GPP_NT_GSL,
        p_eos25_vs_GPP_DT_GSL,p_eos25_vs_NEP_GSL,
        p_cum_eos25_vs_GPP_NT_GSL,
        p_cum_eos25_vs_GPP_DT_GSL,p_cum_eos25_vs_NEP_GSL)
#adding more for more periods:
p_eos90_vs_XEP_GSLsol90<-plot_grid(p_eos90_vs_GPP_NT_GSLsol90,
                                  p_eos90_vs_GPP_DT_GSLsol90,p_eos90_vs_NEP_GSLsol90,
                                  p_cum_eos90_vs_GPP_NT_GSLsol90,
                                  p_cum_eos90_vs_GPP_DT_GSLsol90,p_cum_eos90_vs_NEP_GSLsol90)
p_eos50_vs_XEP_GSL9050<-plot_grid(p_eos50_vs_GPP_NT_GSL9050,
                                  p_eos50_vs_GPP_DT_GSL9050,p_eos50_vs_NEP_GSL9050,
                                  p_cum_eos50_vs_GPP_NT_GSL9050,
                                  p_cum_eos50_vs_GPP_DT_GSL9050,p_cum_eos50_vs_NEP_GSL9050)
p_eos25_vs_XEP_GSL5025<-plot_grid(p_eos25_vs_GPP_NT_GSL5025,
                                  p_eos25_vs_GPP_DT_GSL5025,p_eos25_vs_NEP_GSL5025,
                                  p_cum_eos25_vs_GPP_NT_GSL5025,
                                  p_cum_eos25_vs_GPP_DT_GSL5025,p_cum_eos25_vs_NEP_GSL5025)

#save the plots:
# ggsave(p_eos90_vs_XEP_GSL25sol,
#   filename = paste0("./fig/Results_updated/using_filtered_data/ori_data/further_analysis/ori_GSL25sol_lm_plot.png"))
# ggsave(p_eos90_vs_XEP_GSL2590,
#        filename = paste0("./fig/Results_updated/using_filtered_data/ori_data/further_analysis/ori_GSL2590_lm_plot.png"))
# ggsave(p_eos50_vs_XEP_GSL2550,
#        filename = paste0("./fig/Results_updated/using_filtered_data/ori_data/further_analysis/ori_GSL2550_lm_plot.png"))
# ggsave(p_eos25_vs_XEP_GSL,
#        filename = paste0("./fig/Results_updated/using_filtered_data/ori_data/further_analysis/ori_GSL_lm_plot.png"))
#adding more
ggsave(p_eos90_vs_XEP_GSLsol90,
       filename = paste0("./fig/Results_updated/using_filtered_data/ori_data/further_analysis/ori_GSLsol90_lm_plot.png"))
ggsave(p_eos50_vs_XEP_GSL9050,
       filename = paste0("./fig/Results_updated/using_filtered_data/ori_data/further_analysis/ori_GSL9050_lm_plot.png"))
ggsave(p_eos25_vs_XEP_GSL5025,
       filename = paste0("./fig/Results_updated/using_filtered_data/ori_data/further_analysis/ori_GSL5025_lm_plot.png"))



