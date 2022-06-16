rm=(list=ls())
setwd("/cs/usr/gabriellalevy/3d_hackathon")
getwd()

library(ggplot2)
df_alpha<-read.table("RMSDalphaFold.txt")

df_model<-read.table("RMSDourNet.txt")

df<-merge(df_alpha,df_model, "V1")

ggplot(df, aes(x=V2.y,y=V2.x)) + geom_point(color="navy")+ xlim(0,24) + ylim(0,24)+theme_bw()+
geom_abline(intercept = 0, slope = 1)+xlab("SH3-Net")+ylab("AlphaFold2")+ggtitle("RMSD on domain: SH3-Net vs. AlphaFold2")+
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

ggplot(df, aes(x=V3.y,y=V3.x)) + geom_point(color="navy")+ xlim(0,24) + ylim(0,24)+ theme_bw()+
  geom_abline(intercept = 0, slope = 1)+xlab("SH3-Net")+ylab("AlphaFold2")+ggtitle("RMSD on peptide: SH3-Net vs. AlphaFold2")+
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

