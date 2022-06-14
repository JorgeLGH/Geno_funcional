######EJERCICIO 5######
#Selecciona alg ??un experimento que te parezca interesante de micro-arreglosde expresi ??o
#v??ia la herramienta Geo2R de NCBI. Realiza un an ??alisis de exprei ??on diferencial ent
#re al menos dos grupos experimentales que se reporten en el art ??iculo correspondiente.
#Genera el c ??odigo para realizar todos los analisis como las gr ??aficas y almac??enalo en tu cuenta de Github.
#Comenta el codigo de cade secci??on para que me quede claro que entiendes cada paso par ainferir la red.
#Manda la liga de ese c ??odigo. Genera las tablas en formato csv de genes diferencialmente expresados con las siguientes
#condiciones: logFC de al menos 2 y p-value inferior a 0.1.


# Análisis de expresión diferencial con limma
#   Differential expression analysis with limma
#primero se cargaron las libreríias a utilizar 
library(GEOquery) #se instala de BioConductor
library(limma) #previamente instalada 
library(umap) #esta la instalé a parte en R de CRAN
BiocManager :: install ("GEOquery")
# load series and platform data from GEO
#se deben de cargar la serie del GEO con el que se buscó en GEO2R que comiena con GSE
gset <- getGEO("GSE18388", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] #se establecen los parámetros con ayuda de un if else, sus condicionales. 

fvarLabels(gset) <- make.names(fvarLabels(gset)) #se generan columnas para asignar los datos 

gsms <- "11110000" #grupo al que pertenecen cad auno de los miembros de genes 
sml <- strsplit(gsms, split="")[[1]]

#se lleva a cabo un logFC 
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#se agregan ya genes o experimentos a cada grupo realidado, en este caso, control y muted 
gs <- factor(sml)
groups <- make.names(c("control","muted"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # establecer el modelo lineal de los grupos 

#Aquí se seleccionan los grupos de interés para el contraste de los modelos
#Se calculan los coeficientes con base en esa selección 
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

#modelo estadístico con el modelo de Bayes y se establecen parámetros específicos 
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
#se fijan en una tabla con difernetes características para filtrar los datos 
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

#Aquí se generó una tabla de calidad con base en el control y se realiza un histograma para visuallizar los datos
#el histograma se realiza con base en el p valor de 0.1 para ver los genes diferencialmente expresados.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
#observar genes que están sobre regulados, subregulados y no expresados con base en el p valor 
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.1)

#Diagrama de venn vara observar resultados obtenidos hasta ahora 
vennDiagram(dT, circle.col=palette())

# QQ plot estadístico, se quitan todos aquellos grupos que no cumplen con los parámetros--> T
t.good <- which(!is.na(fit2$F))
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (P valor y logFC)
colnames(fit2) #lista de datos
ct <- 1        # grupo contraste 
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2))) #parámetros para el volcano plot 

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.1) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE18388", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# valor de la distribución d elos grupos 
par(mar=c(4,4,2,1))
title <- paste ("GSE18388", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")


ex <- na.omit(ex) # 
ex <- ex[!duplicated(ex), ]  # 
ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # 
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

#promedio de las varianzas para ver precisión de los pesos de los grupos.
plotSA(fit2, main="Mean variance trend, GSE18388")


######################################################################
###############EJERCICIO 6###############################################################################    
#Elabora un programa en R que, mediante una funci ??on que le pida a un
#usuario un user name y le pida un password que ingrese dos veces y verifque que es el mismo password. 
#La informaci ??on debe guardarse de forma acumulativa en un data.frame para cada nuevo registro. Despu ??es de ingre-
#sar esos datos con otra funci ??on pida esos datos y corrobore que tanto el user name como el password son correctos
#AQUI PUES AGREGA QUE SE VAYA A UN DATA FRAME JAJA ESTÁ A LA MITAD ESTE EJERCICIO
agregarUsuario <- function() {
  username <- readline(prompt="Ingresa el nombre de usuario: ")
  pass <- readline(prompt="Ingresa la contraseña: ")
  pass2 <- readline(prompt="Ingresa la contraseña nuevamente: ")
  
  if (pass == pass2) {
    # nuevo dataframe, es dinamico pq los valores de user y pass corresponden a variables
    print("Bienvenido")
  } else {
    print("Las contraseñas no empatan")
  }
}
#Con esta función se pide que primero se agregue el user name y sus 2 contraseñas
# de esta forma, si las dos son correctas, te arroja un mensaje, de lo contrario, te informa que no coinciden 
agregarUsuario() #se corre la función para ver si funciona y efectivamente, te pide tus datos 



################EJERCICIO 7############################################################################### 
#Utiliza el tutorial de WGCNA para generar una red de co expresi ??on de las que se encuentran en dicho recurso web.
#Comenta el c ??odigo de cada secci ??on para que ne quede claro que entiendes cada paso par ainferir la
#red. Exporta un m ??odulo de la red a Cytoscape o a igraph y mejora su representaci ??on visua


#5. Construction of a weighted gene co-expression network and network modules
# Primero se cargan las librerías que se van a necesitar para este ejercicio
library(WGCNA)  
library(cluster)
options(stringsAsFactors = FALSE); #lo que se ejecutará después no será convertido automáticamente en factores; solo si se le indica
load("Simulated-StandardScreening.RData"); #se cargan los datos previamente guardados en las secciones anteriores 
attach(ModuleEigengeneNetwork1) #aquí se van a copiar los datos nuevos 
#=====================================================

ADJ1=abs(cor(datExpr,use="p"))^6 #se realizó una matriz de adyacencia , con beta=6
k=as.vector(apply(ADJ1,2,sum, na.rm=T)) #función para cuando se tienen pocos genes (<500)
k=softConnectivity(datE=datExpr,power=6) #función cuando hay muchos genes
sizeGrWindow(10,5) #condicionales para el histograma 
par(mfrow=c(1,2))
hist(k) #histograma de vectores rrealizados previamente  (K)
scaleFreePlot(k, main="Check scale free topology\n") #se realiza un histograma y un scalee free 
#para observar conectividades y topología 
#====================================================

datExpr=datExpr[, rank(-k,ties.method="first" )<=3600] #aquí se restringe a solamente mostrar los 
#3600 genes más conectados
#====================================================

dissADJ=1-ADJ1 #para clusterización, se requiere una matriz de disimilaridades debido a que en esto 
#se basa un cluster y así poder agruparlos 
#====================================================

dissTOM=TOMdist(ADJ1) #distancia y sobrelapamiento para señalar disimilaridades
collectGarbage()
#====================================================
#basado en la amtriz de adyacencia, se establecen 3 parametros distintos de clusterización, matriz de adyacencia 
pam4=pam(as.dist(dissADJ), 4)
pam5=pam(as.dist(dissADJ), 5)
pam6=pam(as.dist(dissADJ), 6)
#generación de tablas simuladas de los vectores realizados previamente con los distintos parámetros 
table(pam4$clustering, truemodule)  
table(pam5$clustering, truemodule)
table(pam6$clustering, truemodule)
#===================================================
#nuevamente se toman diferentes parámetros de clusterización , con matriz de distancia 
pamTOM4=pam(as.dist(dissTOM), 4)
pamTOM5=pam(as.dist(dissTOM), 5)
pamTOM6=pam(as.dist(dissTOM), 6)
#se generan tablas para observar los parámetros previamente establecidos 
table(pamTOM4$clustering, truemodule)
table(pamTOM5$clustering, truemodule)
table(pamTOM6$clustering, truemodule)
#==================================================

hierADJ=hclust(as.dist(dissADJ), method="average" ) #cluster basado en adyacencia 
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = data.frame(truemodule), dendroLabels = FALSE, hang = 0.03, 
                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )
#se plotea el cluster ordenado jerárquicamente con base en la matriz de adyacencia y sus disimilaridades
#se observa en la gráfica la cantidad que pertenece a cada grupo
#=================================================
#se marcan o colorean todos aquellos que estén dentro d eun mismo parámetro, con base en el plot ya realizado
colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = data.frame(truemodule, colorStaticADJ),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors") #ploteo
#se observa una linea roja de forma horizontal, lo que indica todos aquellos que están dentro de lo establecido
#estos parámetros los decide cada quien con base en sus necesidades 
#===============================================

branch.number=cutreeDynamic(hierADJ,method="tree") #se utiliza al dendrograma realizado como entrada para esta función
# This function transforms the branch numbers into colors, por método tree
colorDynamicADJ=labels2colors(branch.number )#
#==============================================
#con método híbrido 
colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ, 
                                                  cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))
#Se plotean los resultados obtenidos,sin embargo, se muestran ambos métodos en uno mismo (tree+hibrido)
#híbrido: mayor sensibilidad y poca especificidad; estadístico: alta especificidad y baja sensibilidad
sizeGrWindow(10,5)
plotDendroAndColors(dendro = hierADJ, 
                    colors=data.frame(truemodule, colorStaticADJ, 
                                      colorDynamicADJ, colorDynamicHybridADJ), 
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Gene dendrogram and module colors") #plotear dendrograma con sus características específicas
#============================================
#Nuevamente se toma en cuenta el dendrograma, parámetros distintos establecidos, estadístico
hierTOM = hclust(as.dist(dissTOM),method="average");
colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                                                    deepSplit=2, pamRespectsDendro = FALSE)) #parámetros establecidos 
sizeGrWindow(10,5)
plotDendroAndColors(hierTOM, 
                    colors=data.frame(truemodule, colorStaticTOM, 
                                      colorDynamicTOM, colorDynamicHybridTOM), 
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity") #ploteo del dendrograma 
#se compara con el método híbrido, nuevamente se observan las curvas de identificaicón de cluster
#finalemnte se deben de guardar los análisis realizados 
