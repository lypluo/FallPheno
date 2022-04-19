#----------------------------
#Aim:In this script, I intend to compare the EOS extracted by Lu et al., 2022
# with the EOS extracted by myself-->so that the differences of the analysis are 
# not source from the large differences between phenos extraction

library(dplyr)
library(tidyverse)
#----------------------------
#(1)load XinChen's dataset and the Phenos I extracted:
#---------------
load(file = paste0("./data/data_prep/Add_analysis_based_Lu2022/Data_from_Lu2022.RDA"))
df_Lu2022<-df_merge
#
load(paste0("./data/data_prep/Extracted_Phenos/df_VIs_Phenos_for_Lu2022sites.rda"))
df_this<-Phenos_all

#----------------
#(2)comparsion
#----------------
#---a.
df_this_sel<-df_this %>%
  dplyr::select(sitename,Year,trs_eos90,trs_eos75,trs_eos50,trs_eos25)
df_Lu2022_sel<-df_Lu2022 %>%
  dplyr::select(sitename,Year,EOS_1,EOS_2,EOS_3)

df_merge<-left_join(df_this_sel,df_Lu2022_sel,by=c("sitename","Year"))

#---b.plotting
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
#
names(df_merge)
#
df_merge_sel<-df_merge[,c("trs_eos50","trs_eos25","EOS_1","EOS_2","EOS_3")]
M1<-cor(df_merge_sel,use = "complete.obs")
p.mat.M1<-cor.mtest(df_merge_sel)
corrplot(M1, method="color",type = "upper",order = "hclust",tl.col = "black",
         tl.srt = 45,p.mat = p.mat.M1,sig.level = 0.01,addCoef.col = "black",insig = "blank")
#
lm_fit<-lm(df_merge_sel$trs_eos25~df_merge_sel$EOS_2,data = df_merge_sel)
summary(lm_fit)
plot(df_merge_sel$EOS_2,df_merge_sel$trs_eos25,ylim=c(220,360))
abline(lm_fit,col="red")
# abline(0,1,col="blue",lty=2)
