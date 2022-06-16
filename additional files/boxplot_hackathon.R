rm=(list=ls())
setwd("/cs/usr/gabriellalevy/3d_hackathon")
getwd()

library(ggplot2)
df_alpha<-read.table("RMSDforBoxplotalphaFold.txt", header=FALSE)
df_alpha["Model"]="AlphaFold2"

df_model<-read.table("RMSDforBoxplotourNetNew.txt", header=FALSE)
df_model["Model"]="SH3-Net"

df<-rbind(df_alpha,df_model)

p<-ggplot(data=df, aes(x=V2, y=V1))+
  geom_boxplot(aes(fill=Model)) + xlab("") + ylab("RMSD Å")+theme_bw()+
  ggtitle("RMSD on test set")+
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))
p 

