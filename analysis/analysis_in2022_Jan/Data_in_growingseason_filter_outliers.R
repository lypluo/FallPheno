#####################################
#Aim: tidy the data during the growing season(determined by phenophase from VIs)
#####################################
library(dplyr)
#--------------------
#1.load the data
#--------------------
load(paste0("./data/df_VIs_Phenos.rda")) #phenology data
load(paste0("./data/df_GPP_Meteo_andVIs.rda")) #datasets(GPP,Meteo,and VIs)

#--------------------
#2.select the data could be used for analysis 
#--------------------
#select the data for analysis
df.analysis<-c()
for(i in 1:nrow(Phenos_final)){
  Pheno.temp<-Phenos_final[i,]
  df.temp<-df.final %>%
    filter(sitename==Pheno.temp$sitename & Year==Pheno.temp$Year)
  df.analysis<-rbind(df.analysis,df.temp)
}

#a. check the distribution of the Phenos
par(mfrow=c(4,2))
hist(Phenos_final$UD)
hist(Phenos_final$RD)
hist(Phenos_final$trs_sos25)
hist(Phenos_final$trs_sos50)
hist(Phenos_final$trs_sos75)
hist(Phenos_final$trs_eos75)
hist(Phenos_final$trs_eos50)
hist(Phenos_final$trs_eos25)

#b.check the univariate correlation between phenophases:
library(corrplot)
names(Phenos_final)
Phenos.sel<-Phenos_final[,c("UD","prrD","SD","pop","psrD","RD",
  paste0("trs_sos",c(25,50,75)),paste0("trs_eos",c(75,50,25)))]
M<-cor(Phenos.sel,use = "complete.obs")
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
#
p.mat<-cor.mtest(Phenos.sel)
#refer the code here: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogramvo
png(filename = paste("./fig/Phenos_uni_corrlation.png"),width = 768,height = 768)
corrplot(M, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat,sig.level = 0.01,addCoef.col = "black",insig = "blank")
dev.off()

#--------------------
#3.calculate the GPP, NEE daily mean, Cum GPP/NEE during the growing seaseaon in each site-year
#I selected the growing season between [sos50,eos50] at this step
#--------------------
library(lubridate)
df.sum<-c()
for(i in 1:nrow(Phenos_final)){
  Pheno.temp<-Phenos_final[i,]
  if(is.na(Pheno.temp$trs_eos50) | is.na(Pheno.temp$trs_eos50)) break
  df.temp<-df.analysis %>%
    filter(sitename==Pheno.temp$sitename & Year==Pheno.temp$Year) %>%
    mutate(DoY=lubridate::yday(Date)) %>%
    filter(DoY>=Pheno.temp$trs_sos50 & DoY<=Pheno.temp$trs_eos50)
  df.sum.temp<-df.temp %>%
    summarise(sitename=unique(sitename),Year=unique(Year),sos=min(DoY),eos=max(DoY),
              GPP_NT_mean=mean(GPP_NT_VUT_REF,na.rm=T),   #unit:g C m-2 d-1
              GPP_DT_mean=mean(GPP_DT_VUT_REF,na.rm=T),
              NEE_mean=mean(NEE_VUT_USTAR50,na.rm=T)) %>%
    mutate(NEP_mean=c(-NEE_mean),GSL=c(eos-sos+1)) %>% #growing season length
    mutate(GPP_NT_cum=GPP_NT_mean*GSL,GPP_DT_cum=GPP_DT_mean*GSL,
           NEE_cum=NEE_mean*GSL,NEP_cum=NEP_mean*GSL)
  df.sum<-rbind(df.sum,df.sum.temp)
}
#merge the df.sum and Phenology from Phenos_final
# df.sum$Year<-as.factor(df.sum$Year);df.sum$sitename<-as.factor(df.sum$sitename)
# Phenos_final$Year<-as.factor(Phenos_final$Year);
# Phenos_final$sitename<-as.factor(Phenos_final$sitename)
df.merge<-left_join(df.sum,Phenos_final)
#Remove Duplicate Rows based on multiple variables
df.merge<-distinct(df.merge, sitename,Year,.keep_all = TRUE)

#----
#!!filter the outliers-->remove the site-years has negative GPP_NT_mean or GPP_DT_mean
#----
df.outlier<-df.merge %>%
  filter(GPP_NT_mean<=0 |GPP_DT_mean<=0)
#outliers-->US-WCr: 2010; US-Syv:2008;GL-Zah:2006
df.merge.new<-df.merge %>%
  filter(GPP_NT_mean>0 & GPP_DT_mean>0)  ## also remove when GPP value == NA
#test-->have a look the realtionship between GPP and eos
plot(df.merge.new$GPP_NT_mean,df.merge.new$eos)
plot(df.merge.new$GPP_NT_cum,df.merge.new$eos)
#check the sites have low mean GPP
df.lowGPP_NT<-df.merge.new %>%
  filter(GPP_NT_mean < 3)   #site:GL-ZaH many years; CA-Gro:2014
#
plot(df.merge.new$GPP_DT_mean,df.merge.new$eos)
plot(df.merge.new$GPP_DT_cum,df.merge.new$eos)
#check the sites have low mean GPP
df.lowGPP_DT<-df.merge.new %>%
  filter(GPP_DT_mean < 3)   #site:GL-ZaH many years; CA-Gro:2014
#----
#remove the sites that have very low GPP mean
#----
df.merge.final<-df.merge.new %>%
  filter(GPP_NT_mean > 3 & GPP_DT_mean >3)

#----------------------------
#(4)calculate the anomaly of GPP, NEP..
#----------------------------
#
df.annual.mean<-df.merge.final %>%
  group_by(sitename) %>%
  summarise(GPP_NT_multiY_mean=mean(GPP_NT_mean,na.rm=T),
            GPP_DT_multiY_mean=mean(GPP_DT_mean,na.rm=T),
            NEP_multiY_mean=mean(NEP_mean,na.rm=T),
            GPP_NT_multiY_cum=mean(GPP_NT_cum,na.rm=T),
            GPP_DT_multiY_cum=mean(GPP_DT_cum,na.rm=T),
            NEP_multiY_cum=mean(NEP_cum,na.rm=T),
            sos_mean=mean(sos,na.rm=T),
            eos_mean=mean(eos,na.rm=T),
            pop_mean=mean(pop,na.rm=T),
            GSL_mean=mean(GSL,na.rm=T)
            )
#
site.names<-unique(df.merge.final$sitename)
#
df.new<-c()
for(i in 1:length(site.names)){
  sitename_temp<-site.names[i]
  df.annual.mean.temp<-df.annual.mean[df.annual.mean$sitename==sitename_temp,]
  df.temp<-df.merge.final %>%
    filter(sitename==site.names[i]) %>%
    mutate(A_GPP_NT_mean=GPP_NT_mean - df.annual.mean.temp$GPP_NT_multiY_mean,
           A_GPP_DT_mean=GPP_DT_mean - df.annual.mean.temp$GPP_DT_multiY_mean,
           A_NEP_mean=NEP_mean - df.annual.mean.temp$NEP_multiY_mean,
           A_GPP_NT_cum=GPP_NT_cum - df.annual.mean.temp$GPP_NT_multiY_cum,
           A_GPP_DT_cum=GPP_DT_cum - df.annual.mean.temp$GPP_DT_multiY_cum,
           A_NEP_cum=NEP_cum - df.annual.mean.temp$NEP_multiY_cum,
           A_sos=sos - df.annual.mean.temp$sos_mean,
           A_eos=eos - df.annual.mean.temp$eos_mean,
           A_pop=pop - df.annual.mean.temp$pop_mean,
           A_GSL=GSL - df.annual.mean.temp$GSL_mean
           )
  df.new<-rbind(df.new,df.temp)
}

#----------------------------
#(5) fastly evaluate the realtionship
#----------------------------
library(ggplot2)
library(ggpubr)
#------------------------
##I.using all the original values:
#-----------------------
##1.correlation plot-->sos,eos-->based on the trs_sos50 and trs_eos50
#1a. for phenophase, daily average GPP/NEP, and cum GPP/NEP
df.new_sel1<-df.new[,c("sos","eos","pop",
    "GPP_NT_mean","GPP_DT_mean","NEP_mean",
    "GPP_NT_cum","GPP_DT_cum","NEP_cum",
    "GSL")]
M1<-cor(df.new_sel1,use = "complete.obs")
p.mat.M1<-cor.mtest(df.new_sel1)
#refer the code here: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogramvo
png(filename = paste("./fig/using_filtered_data/uni_corrlation_PhenovsGPP.png"),width = 768,height = 768)
corrplot(M1, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M1,sig.level = 0.01,addCoef.col = "black",insig = "blank")
dev.off()

#1b. for phenophase, anomaly daily average GPP/NEP, and anomaly cum GPP/NEP
df.new_sel2<-df.new[,c(
                       # "sos","eos","pop",
                       "A_sos","A_eos","A_pop",
                       "A_GPP_NT_mean","A_GPP_DT_mean","A_NEP_mean",
                       "A_GPP_NT_cum","A_GPP_DT_cum","A_NEP_cum",
                       "A_GSL")]
M2<-cor(df.new_sel2,use = "complete.obs")
p.mat.M2<-cor.mtest(df.new_sel2)
#refer the code here: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogramvo
png(filename = paste("./fig/using_filtered_data/uni_corrlation_anomaly_PhenovsGPP.png"),width = 768,height = 768)
corrplot(M2, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M2,sig.level = 0.01,addCoef.col = "black",insig = "blank")
dev.off()
##2.regression plot
#(1) mean and cum GPP/NEP vs eos
#mean GPP/NEP vs eos
p_eos_vs_GPP_NT<-ggscatter(df.new,
  x="GPP_NT_mean",y="eos",add = "reg.line")+
    stat_smooth(method = "lm",formula = y ~ x)+
    stat_cor(label.x = 1, label.y = 100) +
    stat_regline_equation(label.x = 1, label.y = 72)
lm_temp<-lm(data=df.new,
            eos ~ GPP_NT_mean)
summary(lm_temp)
p_eos_vs_GPP_DT<-ggscatter(df.new,
                           x="GPP_DT_mean",y="eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 1, label.y = 100) +
  stat_regline_equation(label.x = 1, label.y = 72)
p_eos_vs_NEP<-ggscatter(df.new,
                           x="NEP_mean",y="eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 1, label.y = 100) +
  stat_regline_equation(label.x = 1, label.y = 72)

#cum GPP/NEP vs eos
p_cum_eos_vs_GPP_NT<-ggscatter(df.new,
                           x="GPP_NT_cum",y="eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 1, label.y = 100) +
  stat_regline_equation(label.x = 1, label.y = 72)
lm_temp<-lm(data=df.new,
            eos ~ GPP_NT_cum)
summary(lm_temp)
p_cum_eos_vs_GPP_DT<-ggscatter(df.new,
                           x="GPP_DT_cum",y="eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 1, label.y = 100) +
  stat_regline_equation(label.x = 1, label.y = 72)
p_cum_eos_vs_NEP<-ggscatter(df.new,
                        x="NEP_cum",y="eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = 1, label.y = 100) +
  stat_regline_equation(label.x = 1, label.y = 72)
#merge the plots
library(cowplot)
p_ori<-plot_grid(p_eos_vs_GPP_NT,p_eos_vs_GPP_DT,p_eos_vs_NEP,
          p_cum_eos_vs_GPP_NT,p_cum_eos_vs_GPP_DT,p_cum_eos_vs_NEP)
ggsave(p_ori,filename = paste0("./fig/using_filtered_data/ori_lm_plot.png"))
#(2) Anomaly mean and cum GPP/NEP vs eos
#Anomaly mean GPP/NEP vs eos
p_A_eos_vs_GPP_NT<-ggscatter(df.new,
                           x="A_GPP_NT_mean",y="A_eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = -5, label.y = 80) +
  stat_regline_equation(label.x = -5, label.y = 72)
lm_temp<-lm(data=df.new,
            A_eos ~ A_GPP_NT_mean)
summary(lm_temp)

p_A_eos_vs_GPP_DT<-ggscatter(df.new,
                           x="A_GPP_DT_mean",y="A_eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = -3, label.y = 80) +
  stat_regline_equation(label.x = -3, label.y = 72)

p_A_eos_vs_NEP<-ggscatter(df.new,
                        x="A_NEP_mean",y="A_eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = -1, label.y = 80) +
  stat_regline_equation(label.x = -1, label.y = 72)

#Anomaly cum GPP/NEP vs eos
p_A_cum_eos_vs_GPP_NT<-ggscatter(df.new,
                             x="A_GPP_NT_cum",y="A_eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = -500, label.y = 80) +
  stat_regline_equation(label.x = -500, label.y = 72)
lm_temp<-lm(data=df.new,
            A_eos ~ A_GPP_NT_cum)
summary(lm_temp)

p_A_cum_eos_vs_GPP_DT<-ggscatter(df.new,
                             x="A_GPP_DT_cum",y="A_eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = -200, label.y = 80) +
  stat_regline_equation(label.x = -200, label.y = 72)

p_A_cum_eos_vs_NEP<-ggscatter(df.new,
                          x="A_NEP_cum",y="A_eos",add = "reg.line")+
  stat_smooth(method = "lm",formula = y ~ x)+
  stat_cor(label.x = -100, label.y = 80) +
  stat_regline_equation(label.x = -100, label.y = 72)
#merge the plots
p_anomaly<-plot_grid(p_A_eos_vs_GPP_NT,p_A_eos_vs_GPP_DT,p_A_eos_vs_NEP,
          p_A_cum_eos_vs_GPP_NT,p_A_cum_eos_vs_GPP_DT,p_A_cum_eos_vs_NEP)

ggsave(p_anomaly,filename = paste0("./fig/using_filtered_data/anomaly_lm_plot.png"))

