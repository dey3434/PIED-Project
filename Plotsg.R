# Sharmila Dey
# 1 June 2020

load("./pied_grow_coef.rda")

mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               legend.text=element_text(size=11),legend.title=element_text(size=12),
               legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
               axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
               axis.line.x = element_line(color="black", size = 0.3),
               axis.line.y = element_line(color="black", size = 0.3))

pfun_piedgrow<-function(x,temp,pnorm,tnorm,DIA){
  pied_grow_coef[1,1] + pied_grow_coef[2,1]*pnorm + 
    pied_grow_coef[3,1]*tnorm +pied_grow_coef[4,1]*x + 
    pied_grow_coef[5,1]*temp + pied_grow_coef[6,1]*DIA
}

tfun_piedgrow<-function(x,ppt,pnorm,tnorm,DIA){
  pied_grow_coef[1,1] + pied_grow_coef[2,1]*pnorm + 
    pied_grow_coef[3,1]*tnorm +pied_grow_coef[4,1]*ppt + 
    pied_grow_coef[5,1]*x + pied_grow_coef[6,1]*DIA
}

line_type<-c(NA,"solid","dashed","dotted","dotdash")
quartiles<-c(NA,"1st quartile","2nd quartile","3rd quartile","4th quartile")

qp=quantile(grow$ppt_norm)
pnorms=c(0)
for (i in 2:5){
  
  n=mean(c(qp[i-1],qp[i]))
  pnorms[i]<-n
}

qt=quantile(grow$tmp_norm)
tnorms=c(0)
for (i in 2:5){
  
  n=mean(c(qt[i-1],qt[i]))
  tnorms[i]<-n
}

grow$ppt_norm_cat[grow$ppt_norm<=qp[2]]<-quartiles[2]
grow$ppt_norm_cat[grow$ppt_norm>qp[2] & grow$ppt_norm<=qp[3]]<-quartiles[3]
grow$ppt_norm_cat[grow$ppt_norm>qp[3] & grow$ppt_norm<=qp[4]]<-quartiles[4]
grow$ppt_norm_cat[grow$ppt_norm>qp[4] & grow$ppt_norm<=qp[5]]<-quartiles[5]

grow$tmp_norm_cat[grow$tmp_norm<=qt[2]]<-quartiles[2]
grow$tmp_norm_cat[grow$tmp_norm>qt[2] & grow$tmp_norm<=qt[3]]<-quartiles[3]
grow$tmp_norm_cat[grow$tmp_norm>qt[3] & grow$tmp_norm<=qt[4]]<-quartiles[4]
grow$tmp_norm_cat[grow$tmp_norm>qt[4] & grow$tmp_norm<=qt[5]]<-quartiles[5]


piedgrow_plot <- ggplot(data = grow, aes(x = ppt_yr, y = growth, col = tmp_norm_cat)) + 
  geom_point()+
  stat_function(fun=pfun_piedgrow,args=list(temp=mean(grow$tmp_yr),pnorm=pnorms[2],
                                            tnorm=tnorms[2], DIA = mean(grow$DIA)),
                aes(linetype=quartiles[2]),size=0.75,colour="black")+
  stat_function(fun=pfun_piedgrow,args=list(temp=mean(grow$tmp_yr),pnorm=pnorms[3],
                                            tnorm=tnorms[3], DIA = mean(grow$DIA)),
                aes(linetype=quartiles[3]),size=0.75,colour="black")+
  stat_function(fun=pfun_piedgrow,args=list(temp=mean(grow$tmp_yr),pnorm=pnorms[4],
                                            tnorm=tnorms[4], DIA = mean(grow$DIA)),
                aes(linetype=quartiles[4]),size=0.75,colour="black")+
  stat_function(fun=pfun_piedgrow,args=list(temp=mean(grow$tmp_yr),pnorm=pnorms[5],
                                            tnorm=tnorms[5], DIA = mean(grow$DIA)),
                aes(linetype=quartiles[5]),size=0.75,colour="black")+
  scale_colour_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                      values=c("1st quartile"="#41ae76","2nd quartile"="#238b45",
                               "3rd quartile"="#006d2c","4th quartile"="#00441b"),
                      labels=quartiles[2:5])+
  scale_linetype_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                        values=c("1st quartile"="solid","2nd quartile"="dashed",
                                 "3rd quartile"="dotdash","4th quartile"="dotted"),labels=quartiles[2:5])+
  mytheme+labs(x="Annual temperature",y="Growth")
piedgrow_plot


piedgrow_plot2 <- ggplot(data = grow, aes(x = tmp_yr, y = growth, col = tmp_norm_cat)) + 
  geom_point()+
  stat_function(fun=tfun_piedgrow,args=list(ppt=mean(grow$ppt_yr),pnorm=pnorms[2],
                                            tnorm=tnorms[2],DIA = mean(grow$DIA)),
                aes(linetype=quartiles[2]),size=0.75,colour="black")+
  stat_function(fun=tfun_piedgrow,args=list(ppt=mean(grow$ppt_yr),pnorm=pnorms[3],
                                            tnorm=tnorms[3], DIA = mean(grow$DIA)), 
                aes(linetype=quartiles[3]),size=0.75,colour="black")+
  stat_function(fun=tfun_piedgrow,args=list(ppt=mean(grow$ppt_yr),pnorm=pnorms[4],
                                            tnorm=tnorms[4], DIA = mean(grow$DIA)),
                aes(linetype=quartiles[4]),size=0.75,colour="black")+
  stat_function(fun=tfun_piedgrow,args=list(ppt=mean(grow$ppt_yr),pnorm=pnorms[5],
                                            tnorm=tnorms[5], DIA = mean(grow$DIA)),
                aes(linetype=quartiles[5]),size=0.75,colour="black")+
  scale_colour_manual("TMP Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                      values=c("1st quartile"="#41ae76","2nd quartile"="#238b45",
                               "3rd quartile"="#006d2c","4th quartile"="#00441b"),
                      labels=quartiles[2:5])+
  scale_linetype_manual("TMP Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                        values=c("1st quartile"="solid","2nd quartile"="dashed",
                                 "3rd quartile"="dotdash","4th quartile"="dotted"),labels=quartiles[2:5])+
  mytheme+labs(x="Annual temperature",y="Growth")
piedgrow_plot2
