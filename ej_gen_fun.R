#Ejercicios de genómica

#-----ejercicios de langebio----
library(pvclust)
library(affy)
library(mouse4302cdf)
library(vsn)
library(limma)
library(Biobase)
library(hexbin)

#****Cargar los datos****

setwd("C:/Users/fotgo/OneDrive/Documentos/R/Genomica_funcional/mir155_d/mir155/")#poner el directorio deseado
list.files()#ver archivos que estén en el directorio
pd<-read.table("pdata.txt", header = T, as.is = T) #ver tabla de texto que están nombres de archivos y datos básicos
pd
datos_affy<-ReadAffy(filenames = pd$filename)#leer los datos, le doy los nombres de los archivos basado en el original de texto
pd<-pData(datos_affy)#ponerle los datos al objeto que ya tenía separado con nombres y eso
sampleNames(datos_affy) -> pd$name #agregar nombre a los datos que tengo en sus respectivos objetos
datos_affy#correr este objeto nos muestra la información que hay dentro de los mismos, 6 muestras, 45101 genes, chip de raton 4302

#****Análisis de la calidad de los datos*****

#visualizar los datos, pero no sé bien qué está graficando. En teoría se ven las tendencias de los valores
boxplot(datos_affy, col=rainbow(6))#en teoría se están viendo los valores de expresión en ambas
hist(datos_affy, col=rainbow(6))
windows(2) #está chido, deja ver dos gráficas a la vez en distintas pestañas
image(datos_affy[,3])#intensidad en el microarrgelo como tal, es una visualización previa
heatmap(cor(exprs(datos_affy)), symm=T)#ver correlaciones entre las muestras, esperamos que las wt se correlacionen altamente
corClust <- pvclust(exprs(datos_affy), nboot=1, method.dist="correlation")
plot(corClust) #con el dendograma podemos igual checar si se correlacionan o no, aunque una (ko1) tiene  un comportamniento extraño, 
#en las demás muestras, no es problema
pca <- princomp(exprs(datos_affy))# hace un pca con los datos multivariados para observar tendendicias de agrupación
plot(pca$loadings, main="Principal Component Analysis", col=rainbow(6),  pch=19, cex=2)#plotear 
text(pca$loadings, colnames(exprs(datos_affy)), pos=3, cex=0.8)#agregar los nombres de los samples a cada uno

#****Normalización de los datos****

eset<- rma(datos_affy) #usamos algoritmo rma para normalizar los microarreglos
normalize.AffyBatch.methods() #para ver los métodos disponibles para normalizar los datos de los microarreglos

#****calidad post-normalización****

par(mfrow=c(1,2))    # Una misma ventana con multiples-figuras por "row", con 1 file y 2 columnas
boxplot(datos_affy, col=rainbow(6))
boxplot(data.frame(exprs(eset)), col=rainbow(6)) #ver el boxplot de los datos normalizados
par(mfrow=c(1,1))    # Regresar a una figura por ventana
par(mfrow=c(1,2))
corClust <- pvclust(exprs(datos_affy), nboot=1, method.dist="correlation")
plot(corClust, main="Agrupamiento de muestras antes de normalizar")
corClustAfter <- pvclust(exprs(eset), nboot=1, method.dist="correlation")
plot(corClustAfter, main="Agrupamiento de muestras despues de normalizar") #mejor agrupación con los datos normalizados
par(mfrow=c(1,1))
#Una buena estrategia de normalizado debe controlar la variabilidad a lo largo del rango de valores de expresión. Por ejemplo, 
#no queremos que los genes altamente expresados cambien en promedio más que los genes de expresión media. Podemos revisar esto con la siguiente gráfica:
meanSdPlot(exprs(eset))
boxplot(data.frame(exprs(eset)), col="grey")
lines(exprs(eset)["1428027_at",], lwd=2, type="b", col="red")#vemos la expresión de un nodo, como es el que eliminamos, debe de ser congruente 

#****Guardar los resultados****

write.exprs(eset, file="expr_normalizada.txt")#hacer el archivo con la normalizada
save(eset, file = "eset.Rdata") #guardar los objetos de la sesión

####segunda parte

exprs <- read.table("expr_normalizada.txt", header=TRUE, row.names=1)
boxplot(exprs, col=rainbow(6))#checar que datos estén normalizados con boxplot

#****Calcular expresión diferencial****
#queremos hacer un par de matrices, cada una describe  las categorías, lo experimental y las comparaciones (contraste)
types <- factor(c("KO","KO","KO","WT","WT","WT")) #crear unos niveles, los experimentales y los observados
types
design <- model.matrix(~ 0+types)#describe cuáles muestras son las experimentales y cuáles son las de control
colnames(design) <- levels(types)#solo estoy asiganando el nombre de las columnas respectivamente
design
contMatrix <- makeContrasts(KO-WT, levels=design)#aquí definimos como tal que queremos comparar las pruebas, estamos definiendo las comparaciones entre niveles
contMatrix
#ajustar modelos lineales a nuestros datos de expresión, estimar los contrastes deseados y comprimir las varianzas hacia un valor global recomendados por limma
fit  <- lmFit(exprs,design)
fit2 <- contrasts.fit(fit,contMatrix)
fit2 <- eBayes(fit2)
topTable(fit2, number=20, sort.by="p")#aquí ya tenemos los genes expresados diferencialmente, en este caso la tabla nos acomoda los 20 con una 
#expresión diferencial mayor que otros genes (sub o sobreexpresión); se puede pedir el orden con base en distintos parámetros 

#****Anotando los datos de expresión****
#usaremos paquetes de anotación de los genes para poder hacer la identifiación y todo eso. Son paquetes de anotación de Bioconductor de los microarreglos
library(mouse4302.db)
mouse4302()#para ver los datos que vienen en el paquete, los probes que podemos identificar
probes <- names(fit2$sigma)#asignamos los nombres de las sondas que tenemos en nuestro experimento en un objeto
descriptions <- mget(probes,mouse4302GENENAME)#de estas líneas estamos obteniendo distintas etiquetas que ya vienen dentro de paquete de identificación
#y del cual podemos hacer uso para asignar a nuestros probes
symbols <- mget(probes,mouse4302SYMBOL)
entrezids <- mget(probes,mouse4302ENTREZID)
fit2$genes$EntrezID <- unlist(entrezids)#aquí ya estamos asignando los valores que obtuvimos del paquete de identifiación a nuestros genes observados con probes
fit2$genes$Symbol <- unlist(symbols)#unlist funciona porque cada probe solo tiene un valor de cada variable que se está asginando, ej:ID; de otro modo no
fit2$genes$Description <- unlist(descriptions)
topTable(fit2, number=20, sort.by="p")#aquí volvemos a tener los 20 genes, pero hay más información de mejor entendimiento
volcanoplot(fit2, highlight=10, names=fit2$genes$Symbol)#visualizarlo con un volcano plot
deTable <- topTable(fit2, number=nrow(fit2), lfc=log2(1.5), p.value=0.05)#aquí estamos formando una tabla con los genes que, con base en los criterios de
#significancia que estamos estableciendo, son sginificativos
dim(deTable)
fullTable <- topTable(fit2, number=nrow(fit2), sort.by="logFC")#esta es una tabla con todos los datos, solo que en orden deñ cambio en logFC
dim(fullTable)

#****Para guardar los datos****
write.table(fullTable, file="full_results.txt", row.names=FALSE, sep="\t", quote=FALSE)

#****Ejercicios****
#1.-¿Cuántas sondas predicen como diferencialmente expresadas?
deTable
length(deTable$ID.EntrezID)#nos arroja el número de datos en la columna de ID de cada gen, y da un total de 219
#2.-¿Cuántas decrementan y cuántas aumentan su expresión en el KO?
sum(deTable$logFC>0)#aquí muestra cuántos son de mayor expresión
sum(deTable$logFC<0)#aquí muestra cuántos son de menor expresión
#3.-¿Cuántos genes únicos hay en estas listas?
sum(is.na(deTable$ID.EntrezID))#nos da la suma de NA's en la columna de ID, es decir, aquellos que no se tenía registro y que deben de ser "únicos"






#4.-Genera una figura de volcán manualmente, comparando explicitamente logFC contra -log10(P.Value) e incluyendo todas las sondas

plot(fullTable$logFC, fit2$p.value)

#5.-Colorea de rojo las sondas que encontramos como diferencialmente expresadas.
volcanoplot(fit2, highlight=10, names=fit2$genes$Symbol, col=for(i in 1:length(fit2$genes$EntrezID)){
  if(sum(rownames(fit2)[i]==rownames(deTable))>=1){
    "red"
  }else "black"
})
volcanoplot(fit2, highlight=10, names=fit2$genes$Symbol, col=for(i in 1:length(nrow(deTable))){
  ifelse((rownames(fit2)[i]==rownames(deTable)[i]), "red", "black")})


volcanoplot(fit2, highlight=10, names=fit2$genes$Symbol, col=ifelse((rownames(fit2)==rownames(deTable)), "red", "black"))

#6.-Indica con líneas los cortes que se emplearon, tanto de logFC como de P.Value
volcanoplot(fit2, highlight=219,names=fit2$genes$Symbol)+ abline(v=c(0.05,-0.05), h=(1.5))

            