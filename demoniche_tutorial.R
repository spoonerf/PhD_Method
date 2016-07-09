#install.packages("popbio")
#install.packages("demoniche", repos="http://R-Forge.R-project.org")
library(popbio)
library(demoniche)

data("Hmontana")
str(Hmontana)
Hmontana$Niche_values
Hmontana$matrices



noCC_nodispersal<-demoniche_model(modelname = "Hmontana", Niche=FALSE, Dispersal=FALSE, repetitions = 2, foldername = "noCC_nodispersal")
CC_nodispersal<-demoniche_model(modelname = "Hmontana", Niche=TRUE, Dispersal=FALSE, repetitions = 2, foldername = "CC_nodispersal")


dim(noCC_nodispersal) #number of columns within each matrix within noCC_nodispersal

dimnames(noCC_nodispersal) #names of each column within each matrix within noCC_nodispersal

noCC_nodispersal[,"Meanpop", "Matrix_1"]


barplot(cbind(noCC_nodispersal[40,2,], CC_nodispersal[40,2,]), beside=TRUE, legend.text = Hmontana$list_names_matrices, names.arg = c("no Niche values", "with Niche values"))

#population size for simulations under each scenario matrix scenario of human land use in the last year of the simulation. With climate change there was a lower predicted population size

list.files(path="noCC_nodispersal")

load("noCC_nodispersal/population_results.rda")
str(population_results)

load("noCC_nodispersal/eigen_results.rda")
str(eigen_results)

plot(raster(eigen_results$Reference_matrix$sensitivities))


#Parameterising a model with our own data

args(demoniche_setup)


Hmontana$Orig_Populations


####demographic information

library(popbio)
data(hudvrs)
data(hudsonia)
matrices_mine<-cbind(meanmatrix = as.vector(hudmxdef(hudvrs$mean)), sapply (hudsonia, unlist))

head(matrices_mine)

colnames(matrices_mine) <- c("Reference_matrix", "Matrix_1", "Matrix_2", "Matrix_3", "Matrix_4")

matrix(matrices_mine[,"Reference_matrix"], ncol=6, byrow=FALSE)





