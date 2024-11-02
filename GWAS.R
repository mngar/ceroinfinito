#clase MGyGV
#GWAS

# 1. INSTALAR Y CARGAR LOS PAQUETES NECESARIOS
# 2. SETEAR UNA SEMILLA ASI TODOS TIENEN LOS MISMOS DATOS
# 3. SIMULAR MATRIZ DE DATOS GENOMICOS Y FENOTIPICOS (QTLs y EFECTOS)
# 4. PONERLE NOMBRE A LOS INDIVIDUOS Y MARCADORES
# 5. RANKEAR MARCADORES ASOCIADOS Y EFECTOS: EXPORTAR COMO TABLA
# 6. REALIZAR ANALISIS GWAS CON DISTINTAS METODOLOGIAS Y COMPARAR RESULTADOS CON NUESTROS "QTLs REALES"

# instalacion de paquetes
#PAQUETE DE SIMULACION DE DATOS
install.packages("simulMGF")

#PAQUETES PARA CHECKEAR NORMALIDAD DEL FENOTIPO (GRAFICA Y ANALITICAMENTE)
install.packages("ggplot2")
install.packages("nortest")

#PAQUETE DE ANALISIS DE ASOCIACION
install.packages("remotes")
remotes::install_github("jiabowang/GAPIT3")


#SOLO SI LO ANTERIOR NO FUNCIONA PORQUE NO TIENEN CARGADO RTOOLS:
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#CARGA DE PAQUETES
library(simulMGF)
library(ggplot2)
library(nortest)
library(GAPIT)
#library(GAPIT3)

#setear la semilla (para que los resultados puedan reproducirse exactamente)
set.seed(1234)


#SIMULACION DE DATOS
Nind <- 10000
Nmarkers <- 10000
Nqtl <- 50
Esigma <- .5
Pmean <- 25
Perror <- .25

simulN(Nind, Nmarkers, Nqtl, Esigma, Pmean, Perror)
str(nsimout)

#los QTLs
QTL <- cbind(nsimout$QTN, nsimout$Meffects)
QTL <- cbind(QTL,abs(nsimout$Meffects))
colnames(QTL) <- c("marker", "effect", "effabs")
QTL <- as.data.frame(QTL)
QTL <- QTL[order(-QTL$effabs),]

#matriz de datos genomicos
population <- nsimout$geno
colnames(population) <- c(paste("M", 1:Nmarkers,sep = ""))
rownames(population) <- c(paste("IND", 1:Nind,sep = ""))



#matriz de datos fenotipicos
popfeno <- nsimout$pheno
colnames(popfeno) <- "PHENO"
rownames(popfeno) <- c(paste("IND", 1:Nind,sep = ""))
feno <- data.frame ( IND = rownames(popfeno),
                     Pheno = popfeno)
#armamos la matriz de datos genomicos como la requiere GAPIT3
myGD <- as.data.frame(population)
myGD <- cbind(feno$IND,myGD)
colnames(myGD)[1] <- "IND"
#write.csv(myGD, "myGD.csv")

pop1 = myGD[1:200,]
pop2 = myGD[1:1000,]
pop3 = myGD[1:5000,]
pop4 = myGD




#mapa
map <- data.frame (SNP = colnames(population),
                   Chromosome = c(rep(1,(Nmarkers/5)),rep(2,(Nmarkers/5)),rep(3,(Nmarkers/5)),rep(4,(Nmarkers/5)),rep(5,(Nmarkers/5))),
                   Position = c(1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5))                             #Jiabo Wang: GAPIT requiere que el cromosoma 1 tenga >100 marcadores en cromosoma 1.
)

head(map)


#PRUEBAS DE NORMALIDAD
#CHECKEO VISUAL

# HISTOGRAMA Y CUARVA NORMAL: Consiste en representar los datos mediante un histograma
#y superponer la curva de una distribuci�n normal con la misma media y desviaci�n est�ndar
#que muestran los datos.

ggplot(data = feno, aes(x = PHENO)) +
  geom_histogram(aes(y = ..density.., fill = ..count..)) +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C7C7C") +
  stat_function(fun = dnorm, colour = "firebrick",
                args = list(mean = mean(feno$PHENO),
                            sd = sd(feno$PHENO))) +
  ggtitle("Histograma + curva normal te�rica") +
  theme_bw()

#Grafico de cuantiles teoricos (QQplot)
#Consiste en comparar los cuantiles de la distribucion observada con los cuantiles
#teoricos de una distribucion normal con la misma media y desvio estandar que los datos.
#Cuanto mas se aproximen los datos a una normal, mas alineados estan los puntos entorno a la
#recta.
qqnorm(feno$PHENO, pch = 19, col = "gray50")
qqline(feno$PHENO)


#CHECKEO ANALITICO

#test  Lilliefors: modificaci�n del  Kolmogorov-Smirnov.
#El test Lilliefors asume que la media y varianza son desconocidas, estando especialmente desarrollado
#para contrastar la normalidad. Es la alternativa al test de Shapiro-Wilk cuando el n�mero de
#observaciones es mayor de 50. La funci�n lillie.test() del paquete nortest permite aplicarlo.
lillie.test(x = feno$PHENO)

#D es el estadistico que hay que mirar que sea <0.05


##########
## GWAS ##
##########

detach("package:plotly", unload=TRUE)

#GWAS CON POBLACIONES DE DISTINTOS TAMAnOS
myGAPIT <- GAPIT(
  Y=feno,
  GD=pop1,
  GM=map,
  model=c("GLM"),# choose model
  PCA.total=0,                                          # set total PCAs
  Inter.Plot=TRUE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
  PCA.3d=TRUE,                                          # plot 3d interactive PCA
  file.output=T,
  Geno.View.output=FALSE
)

detach("package:plotly", unload=TRUE)

myGAPIT <- GAPIT(
  Y=feno,
  GD=pop2,
  GM=map,
  model=c("GLM"),# choose model
  PCA.total=0,                                          # set total PCAs
  Inter.Plot=TRUE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
  PCA.3d=TRUE,                                          # plot 3d interactive PCA
  file.output=T,
  Geno.View.output=FALSE
)


############ probablemente no de la memoria para los siguientes en RStudio
detach("package:plotly", unload=TRUE)
myGAPIT <- GAPIT(
  Y=feno,
  GD=pop3,
  GM=map,
  model=c("GLM"),# choose model
  PCA.total=0,                                          # set total PCAs
  Inter.Plot=F,                                      # perform interactive plot
  Multiple_analysis=F,                               # perform multiple analysis
  PCA.3d=F,                                          # plot 3d interactive PCA
  file.output=T,
  Geno.View.output=FALSE
)

detach("package:plotly", unload=TRUE)
myGAPIT <- GAPIT(
  Y=feno,
  GD=pop4,
  GM=map,
  model=c("GLM"),# choose model
  PCA.total=0,                                          # set total PCAs
  Inter.Plot=TRUE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
  PCA.3d=TRUE,                                          # plot 3d interactive PCA
  file.output=T,
  Geno.View.output=FALSE
)


#Comparar marcadores asociados con los QTLs "reales"

QTL$SNP <- paste0("M",QTL$marker)

#marcadores pop1

#GLM
GLM1 <- read.csv("GAPIT.Association.PVE.GLM.PHENO_pop1.csv",header = T)
head(GLM1)
aGLM1 <- GLM1[,c(1,4)]
head(aGLM1)
asoc_GLM1 <- merge(QTL, aGLM1, by = "SNP")
asoc_GLM1

#marcadores pop2

#GLM
GLM2 <- read.csv("GAPIT.Association.PVE.GLM.PHENO_pop2.csv",header = T)
head(GLM2)
aGLM2 <- GLM2[,c(1,4)]
head(aGLM2)
asoc_GLM2 <- merge(QTL, aGLM2, by = "SNP")
asoc_GLM2

#marcadores pop3

#GLM
GLM3 <- read.csv("GAPIT.Association.PVE.GLM.PHENO_pop3.csv",header = T)
head(GLM3)
aGLM3 <- GLM3[,c(1,4)]
head(aGLM3)
asoc_GLM3 <- merge(QTL, aGLM3, by = "SNP")
asoc_GLM3

#marcadores pop4

#GLM
GLM4 <- read.csv("GAPIT.Association.PVE.GLM.PHENO_pop4.csv",header = T)
head(GLM4)
aGLM4 <- GLM4[,c(1,4)]
head(aGLM4)
asoc_GLM4 <- merge(QTL, aGLM4, by = "SNP")
asoc_GLM4


#The statistical power of a hypothesis test is the probability of detecting an effect,
#if there is a true effect present to detect.

#Entonces podemos estimar el poder estad�stico de cada evaluaci�n en nuestros datos
#simulados de la siguiente manera (ya que conocemos los marcadores con asociaci�n real
#al caracter, imposible en datos reales):
#pow <- Marcadores_Asociados_Reales/Nqtl

pow_pop1 <- dim(asoc_GLM1)[1]/Nqtl
pow_pop2 <- dim(asoc_GLM2)[1]/Nqtl
pow_pop3 <- dim(asoc_GLM3)[1]/Nqtl
pow_pop4 <- dim(asoc_GLM4)[1]/Nqtl






# Predecir fenotipos en la poblacion de validación a partir de los distintos modelos:

# ajustar un modelo lineal con los marcadores asociados mediante la función glm
# fit.glm <- glm(feno ~ geno[,c(marcadores asociados)])
# por ejemplo:
# colocando entre los parentesis los nombres de los marcadores asociados o las columnas donde se ubican estos en la matriz de genotipos.

fit.glm1 <- glm(feno$PHENO[1:200] ~ population[1:200,asoc_GLM1$marker])
fit.glm1

fit.glm2 <- glm(feno$PHENO[1:1000] ~ population[1:1000,asoc_GLM2$marker])
fit.glm2

fit.glm3 <- glm(feno$PHENO ~ population[,asoc_GLM3$marker])
fit.glm3

fit.glm4 <- glm(feno$PHENO ~ population[,asoc_GLM4$marker])
fit.glm4

fit.glm1$coefficients



pred.GLM1 = fit.glm1$coefficients[[1]]+fit.glm1$coefficients[[2]]*population[,asoc_GLM1$marker[1]]
pred.GLM1

plot(feno$PHENO, pred.GLM1)

pred.GLM2 = fit.glm2$coefficients[[1]]+fit.glm2$coefficients[[2]]*population[,asoc_GLM2$marker[1]]+
  fit.glm2$coefficients[[3]]*population[,asoc_GLM2$marker[2]]+
  fit.glm2$coefficients[[4]]*population[,asoc_GLM2$marker[3]]+
  fit.glm2$coefficients[[5]]*population[,asoc_GLM2$marker[4]]+
  fit.glm2$coefficients[[6]]*population[,asoc_GLM2$marker[5]]+
  fit.glm2$coefficients[[7]]*population[,asoc_GLM2$marker[6]]+
  fit.glm2$coefficients[[8]]*population[,asoc_GLM2$marker[7]]+
  fit.glm2$coefficients[[9]]*population[,asoc_GLM2$marker[8]]+
  fit.glm2$coefficients[[10]]*population[,asoc_GLM2$marker[9]]+
  fit.glm2$coefficients[[11]]*population[,asoc_GLM2$marker[10]]+
  fit.glm2$coefficients[[12]]*population[,asoc_GLM2$marker[11]]+
  fit.glm2$coefficients[[13]]*population[,asoc_GLM2$marker[12]]+
  fit.glm2$coefficients[[14]]*population[,asoc_GLM2$marker[13]]+
  fit.glm2$coefficients[[15]]*population[,asoc_GLM2$marker[14]]+
  fit.glm2$coefficients[[16]]*population[,asoc_GLM2$marker[15]]+
  fit.glm2$coefficients[[17]]*population[,asoc_GLM2$marker[16]]+
  fit.glm2$coefficients[[18]]*population[,asoc_GLM2$marker[17]]+
  fit.glm2$coefficients[[19]]*population[,asoc_GLM2$marker[18]]+
  fit.glm2$coefficients[[20]]*population[,asoc_GLM2$marker[19]]
  
  

  
  pred.GLM2


  plot(feno$PHENO, pred.GLM2)
  

