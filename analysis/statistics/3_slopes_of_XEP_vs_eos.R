#####################################
#Aim: Further evaluate what's the relationship between GPP/NEP mean or cum vs eos
#####################################
library(dplyr)
#--------------------
#1.load the tidied data 
#--------------------
load(paste0("./data/data_prep/Phenos_and_XEP/preprocessed_phenos_andXEP_updated.RDA"))

#----------------------------
#(2) evaluate the realtionship between gpp/nep vs eos for different periods
#caluculate the slopes between gpp/nep vs eos
#----------------------------
library(ggplot2)
library(ggpubr)
#------------------------
##I.using all the original values:
#-----------------------
df_sel<-df.merge.final 
#-----
#write the function to retrive the parametes in the linear regressions:
#-----
extract_lm<-function(df,x_var,y_var){
  # df<-df.proc
  # x_var<-"XEP_GSL25sol_mean"
  # y_var<-"eos90"
  
  lm_out<-lm(df[,y_var]~df[,x_var],data = df)
  lm_sum<-summary(lm_out)
  ##
  t1<-coef(lm_out)
  t2<-confint(lm_out)
  t3<-lm_sum$coefficients[,4]
  t4<-lm_sum$r.squared
  t5<-cor(df[,y_var],df[,x_var],use = "complete.obs")
  
  #merge the parametes:
  t_merge<-c(as.numeric(t1),as.numeric(t2[2,]),
                      as.numeric(t3[2]),t4,t5)
  names(t_merge)<-c("intercept","slope","slope_low","slope_upper","p-value","Rsqure","R")
  t_merge<-c(round(t_merge[1:4],2),round(t_merge[5],4),round(t_merge[6:7],2))
  return(t_merge)  
}

#----
#writing a function to calculate the slopes between XEP and EOS in diff periods:
#----
slopes_fun<-function(df,df_type,math_type){
  # df<-df_sel
  # df_type<-"GPP_NT"
  # math_type<-"cum"

  df.proc<-df %>%
    select(sitename:GSL_eos50eos25,starts_with(df_type) & ends_with(math_type))
  #
  new.names<-gsub(df_type,"XEP",names(df.proc))
  new.names<-gsub(math_type,"Xmath",new.names)
  names(df.proc)<-new.names
  #conduct the linear regression:
  #1) for the GPP_mean:
  lm_eos90_XEP_Xmath_sos25sol<-extract_lm(df.proc,"XEP_GSL25sol_Xmath","eos90")
  lm_eos90_XEP_Xmath_sos2590<-extract_lm(df.proc,"XEP_GSL2590_Xmath","eos90")
  lm_eos50_XEP_Xmath_sos2550<-extract_lm(df.proc,"XEP_GSL2550_Xmath","eos50")
  lm_eos25_XEP_Xmath_sos2525<-extract_lm(df.proc,"XEP_GSL_Xmath","eos25")
  #
  lm_eos90_XEP_Xmath_sol90<-extract_lm(df.proc,"XEP_GSLsol90_Xmath","eos90")
  lm_eos50_XEP_Xmath_9050<-extract_lm(df.proc,"XEP_GSL9050_Xmath","eos50")
  lm_eos25_XEP_Xmath_5025<-extract_lm(df.proc,"XEP_GSL5025_Xmath","eos25")
  #merge:
  lm_sum<-rbind(lm_eos90_XEP_Xmath_sos25sol,lm_eos90_XEP_Xmath_sos2590,
        lm_eos50_XEP_Xmath_sos2550,lm_eos25_XEP_Xmath_sos2525,
        lm_eos90_XEP_Xmath_sol90,lm_eos50_XEP_Xmath_9050,lm_eos25_XEP_Xmath_5025)
  #
  var_names<-rownames(lm_sum)
  t<-gsub("XEP",df_type,var_names)
  tt<-gsub("Xmath",math_type,t)
  rownames(lm_sum)<-tt
  return(lm_sum)
}
##slopes for differnt data types==>mean:
slopes_GPP_NT_mean<-slopes_fun(df_sel,"GPP_NT","mean")
slopes_GPP_DT_mean<-slopes_fun(df_sel,"GPP_DT","mean")
slopes_NEP_mean<-slopes_fun(df_sel,"NEP","mean")
#
slopes_mean<-as.data.frame(rbind(slopes_GPP_NT_mean,
                                 slopes_GPP_DT_mean,slopes_NEP_mean))
##slopes for differnt data types==>cum:
slopes_GPP_NT_cum<-slopes_fun(df_sel,"GPP_NT","cum")
slopes_GPP_DT_cum<-slopes_fun(df_sel,"GPP_DT","cum")
slopes_NEP_cum<-slopes_fun(df_sel,"NEP","cum")
#
slopes_cum<-as.data.frame(rbind(slopes_GPP_NT_cum,
                                slopes_GPP_DT_cum,slopes_NEP_cum))

##---------------
#(3)plotting
##---------------
plot_slopes<-function(df,math_type){
  # df<-as.data.frame(slopes_mean)
  # math_type<-"XEP-mean"
  # 
  #1)add sig.
  for (i in 1:nrow(df)) {
    if(df$`p-value`[i]>=0.05){
      df$sig.[i]<-"na."
    }
    if(df$`p-value`[i]<0.05 & df$`p-value`[i]>=0.01){
      df$sig.[i]<-"*"
    }
    if(df$`p-value`[i]<0.01 & df$`p-value`[i]>=0.001){
      df$sig.[i]<-"**"
    }
    if( df$`p-value`[i]<0.001){
      df$sig.[i]<-"***"
    }
  }
  #2).add x type and data.type
  df$lm_class<-rep(NA,nrow(df))
  if(math_type=="XEP-mean"){
    df[1:14,]$lm_class<-paste0(substr(rownames(df[1:14,]),4,8),"-",
                               substr(rownames(df[1:14,]),22,30))
    df[15:21,]$lm_class<-paste0(substr(rownames(df[15:21,]),4,8),"-",
                                substr(rownames(df[15:21,]),19,28))
  }
  
  if(math_type=="XEP-cum"){
    df[1:14,]$lm_class<-paste0(substr(rownames(df[1:14,]),4,8),"-",
                               substr(rownames(df[1:14,]),21,30))
    df[15:21,]$lm_class<-paste0(substr(rownames(df[15:21,]),4,8),"-",
                                substr(rownames(df[15:21,]),18,28))
  }
  #set the order the factors:
  df$lm_class<-factor(df$lm_class,
          levels=c("eos90-sos25sol","eos90-sos2590","eos50-sos2550",
                   "eos25-sos2525","eos90-sol90","eos50-9050","eos25-5025"))
  #3) add the data type:
  df$data_type<-c(rep("GPP-NT",7),rep("GPP-DT",7),rep("NEP",7))
  #plotting:
  p_slopes<-ggplot(data = df,aes(x=lm_class,y=slope))+
    geom_point()+
    ylab(paste0(math_type," slope"))+
    xlab("")+
    geom_errorbar(aes(ymin=slope_low,ymax=slope_upper),width=.2)+
    geom_hline(yintercept = 0,col="blue",lty=2)+
    facet_wrap(~data_type)+
    theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))
  if(math_type=="XEP-mean"){
    p_slopes<-p_slopes+geom_text(aes(label=sig.,y=slope_upper+0.5))
  }
  if(math_type=="XEP-cum"){
    p_slopes<-p_slopes+geom_text(aes(label=sig.,y=slope_upper+0.02))
  }
  
  return(p_slopes)    
}
#plots:
slopes_XEP_mean<-plot_slopes(slopes_mean,"XEP-mean")
slopes_XEP_cum<-plot_slopes(slopes_cum,"XEP-cum")

#merge the plots
library(cowplot)
#adding more for more periods:
# p_eos90_vs_XEP_GSLsol90<-plot_grid(p_eos90_vs_GPP_NT_GSLsol90,
#                                   p_eos90_vs_GPP_DT_GSLsol90,p_eos90_vs_NEP_GSLsol90,
#                                   p_cum_eos90_vs_GPP_NT_GSLsol90,
#                                   p_cum_eos90_vs_GPP_DT_GSLsol90,p_cum_eos90_vs_NEP_GSLsol90)

#save the plots:
ggsave(slopes_XEP_mean,width = 10,height=6,
       filename = paste0("./manuscript/fig/Results_updated/using_filtered_data/ori_data/slopes_comparisons/sloeps_XEP_mean.png"))
ggsave(slopes_XEP_cum,width = 10,height=6,
       filename = paste0("./manuscript/fig/Results_updated/using_filtered_data/ori_data/slopes_comparisons/sloeps_XEP_cum.png"))



