install.packages("rgbif")
library(rgbif)
pop<-read.csv("LPI_pops_20160523_edited.csv")
head(pop)
plot(pop$Longitude, pop$Latitude)

pop<-pop[pop$Specific_location ==1,]

pop$Binomial<-as.character(pop$Binomial)

bin_bind<-unique(bin_bind)

library(dismo)

#gbif occurrence poiont counting

gbif_count<-function(binomial){
  
  bin<-strsplit(binomial, "_")
  genus<-bin[[1]][1]
  species<-bin[[1]][2]
  sp<-gbif(genus, species, geo=TRUE, download=FALSE,ntries=10)
  records<-cbind(genus, species, sp)
  return(records)
}


gbif_records<-lapply(pop$Binomial, gbif_count)

gbif<-data.frame(matrix(unlist(gbif_records), ncol=3, byrow=T))
gbif$binomial<-paste(gbif$X1, gbif$X2, sep="_")
colnames(gbif)<-c("Genus", "species", "GBIF_count", "Binomial")

gbif_count<-gbif[,-c(1,2)]

gbif_count$GBIF_count<-as.numeric(as.character(gbif_count$GBIF_count))

filter(gbif_count, GBIF_count > 1000000) %>%
  arrange(desc(GBIF_count))

gbif_count<-unique(gbif_count)

#comadre matrix counting
load(paste(wd, "COMADRE_v.2.0.1.RData", sep="/"))


matrix_count<-function(binomial){
  
  species<-binomial
  tempMetadata<-subset(comadre$metadata, SpeciesAccepted==species)
  mat_count<-nrow(tempMetadata)
  mat_sp<-cbind(species, mat_count)
  return(mat_sp)
}


com_count<-lapply(lpi$binomial, matrix_count)
com_count<-data.frame(matrix(unlist(com_count), ncol= 2, byrow=T))

com_count$X2<-as.numeric(as.character(com_count$X2))
colnames(com_count)<-c("Binomial", "Matrix_count")


#lpi pop counting
pop2<-pop %>% group_by(Binomial) %>% mutate(count = n())
pop2$count

lpi_pop<-data.frame(pop2$Binomial, pop2$count, pop2$Red_list_category, pop2$Common_name, pop2$Class, pop2$System)
colnames(lpi_pop)<-c("Binomial", "pop_count", "RedList", "Common_Name", "Class", "System")

gbif_count

com_count


lpi_gbif<-merge(lpi_pop, gbif_count, by="Binomial")

lgc<-merge(lpi_gbif, com_count, by="Binomial")

lgc<-unique(lgc)

head(lgc)

lgc$Matrix_count<-as.numeric(as.character(lgc$Matrix_count))

lgc_ord<-lgc %>%
  filter(Matrix_count >1 & pop_count > 5 & GBIF_count>1000 & System == "Terrestrial") %>%
  arrange(desc(Matrix_count), desc(pop_count))

lgc_ord

write.csv(lgc_ord, "Comadre_gbif_lpi_counts.csv")





