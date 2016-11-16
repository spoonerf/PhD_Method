bm<-read.csv("bodymass_missing.csv")
birds<-read.csv("Birdlife_Bodymass.csv")
mamms<-read.csv("Pantheria_update_8.6.05.csv")

birds<-data.frame(birds$Binomial_BL, birds)


no_bm_b<-subset(bm,is.na(Bodymass_g)&Class=="Aves")

no_bm_m<-subset(bm,is.na(Bodymass_g)&Class=="Mammalia")

no_bm_b$Binomial %in% birds$Binomial_BL

bm_b<-merge(no_bm_b, birds, by.x="Binomial", by.y ="Binomial_BL")

subset(bm_b, is.na(Bodymass.x))
