all_year_ab<-read.csv("all_year_ab_ursus.csv")
all_year_ab$ID<-as.numeric(as.character(all_year_ab$ID))

sdm_lpi_melt<-read.csv("sdm_melt_ursus.csv")

melt_lambda_short<-read.csv("melt_lambda_short_ursus.csv")


melt_lambda_short$run_id<-paste(melt_lambda_short$ID, melt_lambda_short$ldd, melt_lambda_short$SD, melt_lambda_short$rep_id)

sp_lam_avg<-melt_lambda_short %>%
  group_by(ID, Year)%>%
  summarise(mean_lambdas = mean(Lambdas))



melt_lambda_short<-filter(melt_lambda_short, Year >= 1970 & ID != 4357 & ID != 4358)
sp_lam_avg<-filter(sp_lam_avg, Year >= 1970 & ID != 4357 & ID != 4358)
all_year_ab<-filter(all_year_ab, Year >= 1970 & ID != 4357 & ID != 4358)
sdm_lpi_melt<-filter(sdm_lpi_melt, Year >= 1970 & ID != 4357 & ID != 4358)

sdm_lpi_melt<-sdm_lpi_melt %>%
  group_by(ID)%>%
  mutate(log_hsi = log10(HSI))%>%
  mutate(HSI_Lambdas = c(diff(log_hsi), NA))

ggplot()+
  geom_line(data = melt_lambda_short, aes(x= Year, y=Lambdas, group=run_id), colour = "grey")+
  geom_line(data = sp_lam_avg, aes(x = Year, y = mean_lambdas), colour = "black", size = 1)+
  geom_line(data = sdm_lpi_melt, aes(x = Year, y = HSI_Lambdas), colour = "blue", size = 1)+
  geom_line(data = all_year_ab, aes(x = Year, y = Lambdas), colour = "red", size = 1.5)+
  theme_bw()+
  facet_wrap(.~ ID)



###########Red deer##############

all_year_ab<-read.csv("all_year_ab_cervus.csv")
all_year_ab$ID<-as.numeric(as.character(all_year_ab$ID))

sdm_lpi_melt<-read.csv("sdm_melt_cervus.csv")

melt_lambda_short<-read.csv("melt_lambda_short_cervus.csv")

melt_lambda_short$run_id<-paste(melt_lambda_short$ID, melt_lambda_short$ldd, melt_lambda_short$SD, melt_lambda_short$rep_id)

sp_lam_avg<-melt_lambda_short %>%
  group_by(ID, Year)%>%
  summarise(mean_lambdas = mean(Lambdas))



melt_lambda_short<-filter(melt_lambda_short, Year >= 1970 & ID == 542 | ID == 6555 | ID == 11176 | ID == 11180)
sp_lam_avg<-filter(sp_lam_avg, Year >= 1970 & ID == 542 | ID == 6555 | ID == 11176 | ID == 11180)
all_year_ab<-filter(all_year_ab, Year >= 1970 & ID == 542 | ID == 6555 | ID == 11176 | ID == 11180)
sdm_lpi_melt<-filter(sdm_lpi_melt, Year >= 1970 & ID == 542 | ID == 6555 | ID == 11176 | ID == 11180)

sdm_lpi_melt<-sdm_lpi_melt %>%
  group_by(ID)%>%
  mutate(log_hsi = log10(HSI))%>%
  mutate(HSI_Lambdas = c(diff(log_hsi), NA))

ggplot()+
  geom_line(data = melt_lambda_short, aes(x= Year, y=Lambdas, group=run_id), colour = "grey")+
  geom_line(data = sp_lam_avg, aes(x = Year, y = mean_lambdas), colour = "black", size = 1)+
  geom_line(data = sdm_lpi_melt, aes(x = Year, y = HSI_Lambdas), colour = "blue", size = 1)+
  geom_line(data = all_year_ab, aes(x = Year, y = Lambdas), colour = "red", size = 1.5)+
  theme_bw()+
  facet_wrap(.~ ID)



##########Alpine Ibex


all_year_ab<-read.csv("all_year_ab_capra_new.csv")
all_year_ab$ID<-as.numeric(as.character(all_year_ab$ID))

sdm_lpi_melt<-read.csv("sdm_melt_capra.csv")

melt_lambda_short<-read.csv("melt_lambda_short_capra_new.csv")

melt_lambda_short$run_id<-paste(melt_lambda_short$ID, melt_lambda_short$ldd, melt_lambda_short$SD, melt_lambda_short$rep_id)

sp_lam_avg<-melt_lambda_short %>%
  group_by(ID, Year)%>%
  summarise(mean_lambdas = mean(Lambdas))



melt_lambda_short<-filter(melt_lambda_short, Year >= 1970 )
sp_lam_avg<-filter(sp_lam_avg, Year >= 1970)
all_year_ab<-filter(all_year_ab, Year >= 1970)
sdm_lpi_melt<-filter(sdm_lpi_melt, Year >= 1970)

sdm_lpi_melt<-sdm_lpi_melt %>%
  group_by(ID)%>%
  mutate(log_hsi = log10(HSI))%>%
  mutate(HSI_Lambdas = c(diff(log_hsi), NA))

ggplot()+
  geom_line(data = melt_lambda_short, aes(x= Year, y=Lambdas, group=run_id), colour = "grey")+
  geom_line(data = sp_lam_avg, aes(x = Year, y = mean_lambdas), colour = "black", size = 1)+
  geom_line(data = sdm_lpi_melt, aes(x = Year, y = HSI_Lambdas), colour = "blue", size = 1)+
  geom_line(data = all_year_ab, aes(x = Year, y = Lambdas), colour = "red", size = 1.5)+
  theme_bw()+
  facet_wrap(.~ ID)


