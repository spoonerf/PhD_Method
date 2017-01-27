betts_2005<-0.76
ipcc_ar5_1986_2005<-0

model_rate_birds<-0.0701
model_pop_birds<- -0.0474
model_pop_birds_lci<- -0.0251
model_pop_birds_uci<- -0.0692  

model_rate_mammals<-0.079
model_pop_mammals<- -0.0207
model_pop_mammals_lci<- -0.0066
model_pop_mammals_uci<- -0.0346  

pop_extrap<-function(year,temp, class){
  year_diff<-year - 2005
  temp_diff<-temp - ipcc_ar5_1986_2005#betts_2005
  rate<-temp_diff/year_diff
  
  if (class=="Birds"){
   rate_birds<-rate/model_rate_birds  
   lambda_mean<-rate_birds*model_pop_birds
   lambda_mean_lci<-rate_birds*model_pop_birds_lci
   lambda_mean_uci<-rate_birds*model_pop_birds_uci
   
   pop_decline<-(1 + lambda_mean)^year_diff
   pop_decline_lci<-(1 + lambda_mean_lci)^year_diff
   pop_decline_uci<-(1 + lambda_mean_uci)^year_diff
  
   pop_decline_int<-c(pop_decline_lci, pop_decline, pop_decline_uci) 
    } else {
      rate_mammals<-rate/model_rate_mammals  
      lambda_mean<-rate_mammals*model_pop_mammals
      lambda_mean_lci<-rate_mammals*model_pop_mammals_lci
      lambda_mean_uci<-rate_mammals*model_pop_mammals_uci
      
      pop_decline<- (1 + lambda_mean)^year_diff
      pop_decline_lci<-(1 + lambda_mean_lci)^year_diff
      pop_decline_uci<-(1 + lambda_mean_uci)^year_diff
      
      pop_decline_int<-c(pop_decline_lci, pop_decline, pop_decline_uci) 
      }

  return(pop_decline_int)
}

# years<-seq(2040,2100,10)
# temps<-rep(1.5, 7)
# birds<-rep("Birds",7)


RCP2.6_blci<-pop_extrap(2100,0.3, "Birds")
RCP2.6_buci<-pop_extrap(2100,1.7, "Birds")

RCP4.5_blci<-pop_extrap(2100,1.1, "Birds")
RCP4.5_buci<-pop_extrap(2100,2.6, "Birds")

RCP6.0_blci<-pop_extrap(2100,1.4, "Birds")
RCP6.0_buci<-pop_extrap(2100,3.1, "Birds")

RCP8.5_blci<-pop_extrap(2100,2.6, "Birds")
RCP8.5_buci<-pop_extrap(2100,4.8, "Birds")


RCP2.6_mlci<-pop_extrap(2100,0.3, "Mammals")
RCP2.6_muci<-pop_extrap(2100,1.7, "Mammals")

RCP4.5_mlci<-pop_extrap(2100,1.1, "Mammals")
RCP4.5_muci<-pop_extrap(2100,2.6, "Mammals")

RCP6.0_mlci<-pop_extrap(2100,1.4, "Mammals")
RCP6.0_muci<-pop_extrap(2100,3.1, "Mammals")

RCP8.5_mlci<-pop_extrap(2100,2.6, "Mammals")
RCP8.5_muci<-pop_extrap(2100,4.8, "Mammals")

Class<-rep(c("Bird", "Mammal"), each=4)
rcp<-rep(c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"),2)
lci<-c(RCP2.6_blci[2],RCP4.5_blci[2],RCP6.0_blci[2],RCP8.5_blci[2],RCP2.6_mlci[2],RCP4.5_mlci[2],RCP6.0_mlci[2],RCP8.5_mlci[2] )
uci<-c(RCP2.6_buci[2],RCP4.5_buci[2],RCP6.0_buci[2],RCP8.5_buci[2],RCP2.6_muci[2],RCP4.5_muci[2],RCP6.0_muci[2],RCP8.5_muci[2] )


df<-data.frame(Class,rcp,lci,uci)
df$lci<-df$lci*100
df$uci<-df$uci*100

library(ggplot2)
q1<-ggplot(df, aes(colour=Class))
q1<- q1 + geom_linerange(aes(x=rcp, ymin=lci, ymax=uci), lwd=2.5, position = position_dodge(width=2/3))
q1<- q1 + theme_bw() + labs(y = "Remaining Population in 2100 (%)", x = "") + theme(legend.position="none",text=element_text(size=20),axis.text.x=element_text(size=20) , axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
q1<- q1 + scale_y_continuous(limits=(c(0,100)))
q1<- q1 + scale_color_manual(values=c("black", "black"))
q1







pop_extrap(2100,1.5,"Birds")
pop_extrap(2100,2.0,"Birds")
pop_extrap(2100,4.0,"Birds")
pop_extrap(2100,6.0,"Birds")

pop_extrap(2050,1.5,"Birds")
pop_extrap(2050,2.0,"Birds")
pop_extrap(2050,4.0,"Birds")
pop_extrap(2050,6.0,"Birds")


pop_extrap(2100,1.5,"Mammals")
pop_extrap(2100,2.0,"Mammals")
pop_extrap(2100,4.0,"Mammals")
pop_extrap(2100,6.0,"Mammals")

pop_extrap(2050,1.5,"Mammals")
pop_extrap(2050,2.0,"Mammals")
pop_extrap(2050,4.0,"Mammals")
pop_extrap(2050,6.0,"Mammals")


