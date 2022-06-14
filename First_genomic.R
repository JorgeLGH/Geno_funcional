library(ggplot2)
library(Biostrings)
library(BSgenome)
library(igraph)
#notes:
#Biological networks:
#Introduction:
#impacto de redes a nivel mundial, no solo en ciencias biológicas
#últmios años ha sido el boom
#utilidad: se aplican en muchas áreas, no solo en ciencias o ciencias biológicas
#cosas que están conectadas forman una red
#principalmente impulsado por el ejército por ejército para evitar ataqiues terroristas

###Área definida de ciencia, COMPLEJIDAD
#Se refiere a algo con muchas partes interconectadas, o puede ser algo complejo
#las redes nos ayudan a entender la complejidad
#es una ciencia emergente y con muchas aplicaciones
##es una teoría científica que propone que algunos sistemas se comportan en un modo tal que,
#no es posible entenderlo a través de análsis estándar
#fenómenos emergentes son mejor entendidos a través de interacciones de compuestos internos
#de un sistema
##Propiedad emergente:
#existen comportamientos que emergen a partir de la combinación del conjunto de varios
#sistemas interactuandoen un ambiente, formando comportamientos más complejos en colectivo
#Un sistema es complejo es más que la suma de sus partes



####Exercise figures#####

####1.-10 nodos
##Degree
#nodo 1---- 9
#resto----- 1
##Av. degree
#1.8
##Degree distribution
#0.9---1 connection
#0.1---9 connections

#2.-Ring
#Degree
#todos 2
#Av. degree
#2
#Degree distribution
#100%---2

#3.-Above Right
#Degree
#1,2,5----5
#3------3
#4,7,8,9------6
#6,10------4
#av. degree
#5
#degree distribution
#0.3----5
#0.1----3
#0.4----6
#0.2----4

#4.-down right
#Degree
#1------5
#2,3------2
#4,6,7,8,9,10------1
#5------3
#Av. degree
#1.8
#Degree distribution
#0.1----5
#0.2----2
#0.6----1
#0.1----3


###---igrpah_first_exercises----

g_1<-make_empty_graph(n=10, directed = T)#make an empty network
g_1
V(g_1)$color="blue"#color of the shape
V(g_1)$shape="square"#shape

plot(g_1)
g_1<-add.edges(g_1, c(1,3, 1,4, 1,7, 2,3, 2,5, 3,9, 4,8, 6,2, 6,10))#add the connections
g_1<-add.vertices(g_1, 2, color="green", shape="square")#may add another "point"
plot(g_1)

class(g_1)
str(g_1)

g_2<-delete.edges(g_1, c(1))#delete the connection with the selected number
plot(g_2)
g_2<-add.edges(g_2, c(1,10, 11,12))#add any edge you want
plot(g_2)


V(g_2)$name<- letters[10:21]
E(g_1)#shows the connections
V(g_2)#shows the names of the nodes and # of vertex

#Degree
degree(g_2)#shows the degree of each node
plot(g_2, layout=layout_nicely, vertex.size=degree(g_2, V(g_2), "in")*15+20,
     vertex.label.dist=0.5, edge.arrow.size=.5)
hist(degree(g_2), col = 4)

#Degree distribution
plot(degree.distribution(g_2), main="Degree Distribution", 
     xlab="Degree", ylab="Frequency", col=c(2:5), pch=20, lwd=4)


#Matrix
adj_mat<-as.matrix(get.adjacency(g_2))
adj_mat
get.adjacency(g_2)

#Default (may use to practice)
g_3<-make_star(20)
g_3
plot(g_3)
plot(make_ring(10))
g_4<-make_ring(10, directed = T)
g_4
plot(g_4)
plot(make_kautz_graph(2,3))

#heatmap
heatmap(adj_mat, Rowv=NA, Colv="Rowv")

##----Exercises_with_classmates----
library(igraph)
fr<-read.csv("RedDeAmistadesGF2021Nombres - Hoja 1.csv")
rownames(fr)<-fr[,1]
fr_mat<-as.matrix(fr[,-1])
diag(fr_mat)<-rep(0,19)
red_1<-graph_from_adjacency_matrix(fr_mat, mode= "directed")
red_1
###
plot(red_1)
degree(red_1, mode="in")
plot(cluster_spinglass(red_1),red_1)
hist(degree(red_1,mode="in"), col=c(1:length(hist(degree(red_1, mode="in")))))
hist(degree(red_1,mode="out"), col=c(1:length(hist(degree(red_1, mode="out")))))

##----Second_ex_friend_network----

#Distances

shortest.paths(red_1)#shortest path between a node and the others
distances(red_1, v= "ANA", to="ENOE")#may select a single node to other
diameter(red_1)
mean_distance(red_1)

plot(red_1,layout=layout.fruchterman.reingold(red_1), 
     vertex.size=5, edge.arrow.size=.2)
###
g <- barabasi.game(1000, power=1)
plot(g,layout=layout.fruchterman.reingold(g), 
     vertex.size=5, edge.arrow.size=.2)
###

mat.1<-shortest.paths(red_1)#apparently seems to be a matrix by default
heatmap(mat.1)
shortest_paths(red_1,1,6)
all_shortest_paths(red_1, 1,2:6)




#-------exercise_find_solution-----
for (i in 1:length(a)) {
  if(sum(a[i]>a[1:length(a)])==(length(a)-1)){
    print(paste(a[i], "es el más grande"))
  }
}
a<-c(13,21,54,42,3,54)
##
for (i in 1:length(a)) {
  if(sum(a[i]>=a[1:length(a)])==(length(a))){
    print(paste(a[i], "es el más grande"))
  }
}
a<-c(13,21,54,42,3,54,100)
  
#----idk-----
alt_h<-c(1.70, 1.73, 1.82, 1.76, 1.75, 1.76)
alt_m<-c(1.63, 1.59, 1.65, 1.69, 1.70, 1.60, 1.68, 1.64, 1.59, 1.58, 1.63)
hist(alt_h)
hist(alt_m)
sd(alt_h)
sd(alt_m)
mean(alt_h)
mean(alt_m)

red_1

barplot(sort(degree(red_1, mode = "in")), beside = T, col = rainbow(19))
barplot(sort(degree(red_1, mode = "out")), beside = T, col = rainbow(19))

plot(sort(degree(red_1, mode = "in")))
plot(sort(degree(red_1, mode = "out")))

sort(degree(red_1, mode = "in"), decreasing = T)[1:3]
sort(degree(red_1, mode = "out"), decreasing = T)[1:3]
####

red_1["JORGE"]
as.vector(red_1["JORGE"])
aa<-as.vector(names(which(red_1["JORGE"]==1)))
ma.1<-as.matrix(red_1[aa])
mean(rowSums(ma.1))
#

#------robustness

g100<-barabasi.game(100,directed = FALSE)
g1K<-barabasi.game(1000,directed = FALSE)
g2<-random.graph.game(1000,0.20)
g3<-sample_smallworld(1,1000,p=0.2,nei=3)
##
rem.rd.1<-function(z){
  for (i in 1:10) {
    z<-delete.vertices(z,sample(1:length(V(z)),1))
    print(paste("the mean distance for",i, "deletions is", mean_distance(z)))
    print(paste("the diameter for",i, "deletions is", diameter(z)))
  }
}

rem.rd.1(g100)
rem.rd.1(g1K)
rem.rd.1(g2)
rem.rd.1(g3)

##
rem.sl.1<-function(z){
  for (i in 1:10) {
    f<-delete.vertices(z,sort(degree(z),decreasing = T)[1])
    print(paste("the mean distance for",i, "deletions is", mean_distance(f)))
    print(paste("the diameter for",i, "deletions is", diameter(f)))
  }
}

rem.sl.1(g100)
rem.sl.1(g1K)
rem.sl.1(g2)
rem.sl.1(g3)


#-----

library(igraph)

a<-barabasi.game(100, directed = F)
degree(a)
sort(degree(a), decreasing = T)
eccentricity(a)
sort(eccentricity(a),decreasing = F)
closeness(a)
betweenness(a)


sort(degree(a), decreasing = T)
sort(eccentricity(a),decreasing = F)
sort(closeness(a),decreasing = T)
sort(betweenness(a),decreasing = T)
eigen_centrality(a)
page.rank(a)
#with my network
sort(degree(red_1, mode = "in"), decreasing = T)
sort(degree(red_1, mode = "out"), decreasing = T)
sort(eccentricity(red_1, mode = "in"),decreasing = F)
sort(eccentricity(red_1, mode = "out"),decreasing = F)
sort(closeness(red_1, mode = "in"),decreasing = T)
sort(closeness(red_1, mode = "out"),decreasing = T)
sort(betweenness(red_1, directed = T),decreasing = T)
eigen_centrality(red_1, directed = T)#creo que es en relación de los nodos con más conexiones, esos nodos qué tan conectados con los demás
page.rank(red_1)


plot

library(igraphdata)
####HW

#-----
mm<-matrix(c(0,1,1,1,0,0,1,0,0), ncol = 3, nrow = 3)
rr<-graph_from_adjacency_matrix(mm, mode = "undirected")
V(rr)$name<-letters[1:3]
eigen_centrality(rr)[1]
eigen(mm)
#----
library(igraphdata)
data(package="igraphdata")
data(yeast)
yeast
#----
aa<-edge.betweenness.community(red_1,directed = T, modularity = NULL)
bb<-label.propagation.community(red_1, weights = NULL)
cc<-cluster_optimal(red_1, weights = NULL)
dd<-cluster_spinglass(red_1)
plot(aa, red_1)
plot(bb, red_1)
plot(cc, red_1)
plot(dd, red_1)
plot(leading.eigenvector.community(red_1),red_1)
#-----
log(10000100000)
log(1000000)
log(log(10000100000))
log(log(1000000))
aa<-barabasi.game(100, directed = F)
bb<-barabasi.game(1000, directed = F)
cc<-random.graph.game(100, 1/10)
dd<-random.graph.game(1000, 1/10)
ee<-sample_smallworld(1, 100, 5, 0.1, 0.05)
diameter(aa)
diameter(bb)
mean(distances(aa))
mean(distances(bb))
diameter(cc)
mean(distances(cc))
diameter(dd)
mean(distances(dd))
diameter(ee)
mean(distances(ee))

library(igraph)
library(get.adj)
kk<-read.csv("RedDeAmistadesGF2021Nombres_nueva")
kk
kk_mat<-as.matrix(kk)
diag(kk_mat)<-rep(0,15)
red_1<-graph_from_adjacency_matrix(kk_mat, mode= "directed")
red_1
get.adjlist(red_1, mode = "out")
ff<-as.matrix(as_edgelist(red_1, names = T))
write.csv(ff, "prueba.csv")

#-----
#redes de gustos 
gg<-read.csv("Red de co-gustos - Hoja 1.csv")
rownames(gg)<-gg[,1]
gg<-gg[,-1]
gg<-t(as.matrix(gg))#transponer
cor(gg, method = "pearson")
cor(gg, method = "kendall")
cor(gg, method = "spearman")
library(pheatmap)
pheatmap(cor(gg, method = "pearson"))
#establecer criterio de corte para la red, de modo que sea para ver si se conecta el nodo o no
#necesitamos cambiar los valores negativos, se le suma 1 para que todo de positivo
mm<-cor(gg, method = "pearson")
mm<-(mm+1)/2
#valor arbitrario, en este caso 0.7
diag(mm)<-rep(0,15)
for (i in 1:length(mm)) {
  if(mm[i]>0.7){
    mm[i]<-1
  }else{
    mm[i]<-0
  }
}
ff<-graph.adjacency(mm,mode = "undirected")
write.csv(mm, "gustos_1.csv")
rr<-graph_from_adjacency_matrix(mm, mode="undirected", weighted=T)
plot(cluster_spinglass(rr, weights = E(rr)$weight),rr)
##Formar red pesada
library(WGCNA)
library(igraphinshiny)
library(igraph)
rr
plot(rr,edge.width=E(rr)$weight*1, edge.color="blue")
write.csv(mm, "pesada_1.csv")
hh<-matrix(as_edgelist(rr, names = T), ncol = 2)


#-----
library(BoolNet)
red_3<-loadNetwork("Primera_red.txt")
plotNetworkWiring(red_3)
atrac<-getAttractors(red_3)
plotAttractors(atrac)
plotSequence(red_3, includeAttractorStates = "all", startState = rep(1,2))#not working
getPathToAttractor(red_3, state = c(0,1,1), includeAttractorStates = "all")

red_pina<-loadNetwork("new_2.txt")
plotNetworkWiring(red_pina)
atrac<-getAttractors(red_3)
plotAttractors(atrac)

#------
library(igraph)
library(WGCNA)
clincal_traits<-read.csv("FemaleLiver-Data/ClinicalTraits.csv")
gene_annotation<-read.csv("FemaleLiver-Data/GeneAnnotation.csv")
liver_female<-read.csv("FemaleLiver-Data/LiverFemale3600.csv")
#limpieza de los datos
options(stringsAsFactors = FALSE)
femData<-liver_female
dim(femData)
names(femData)
datExpr0 <- as.data.frame(t(femData[, -c(1:8)]))
names(datExpr0) <- femData$substanceBXH
rownames(datExpr0) <- names(femData)[-c(1:8)]

#----

library(BoolNet)
data(cellcycle)
plotStateGraph(getAttractors(cellcycle))

#-----
library(igraph)
library(BoolNet)
lac<-loadNetwork("lac_bool.txt")
lac
plotNetworkWiring(lac)
ff<-generateState(lac, c("M"=1, "P"=1))
getAttractors(lac, startStates = list(ff))

#-----ejercicios de langebio----
library(pvclust)
library(affy)
library(mouse4302cdf)
library(vsn)
library(limma)
library(Biobase)
library(hexbin)

#****Cargar los datos****

setwd("mir155_d/mir155/")#poner el directorio deseado
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
design <- model.matrix(~ 0+types)
colnames(design) = levels(types)
design








rownames(fit2[[1]])->probes
probes <- names(fit2$sigma)
#repositorio de archivos cel

geo



#----
BiocManager::install("rhdf5")
devtools::install_github("pachterlab/sleuth")
library(devtools)
library(rhdf5)
library(sleuth)
library(edgeR)
library(ensembldb)

setwd("C:/Users/fotgo/OneDrive/Documentos/R/Genomica_funcional/")

#----coso----

#10 soles, 5 y 5, 3 y 7, 1 y 9

binomial.prob<-function(x,y){
  probs<-readline(prompt=paste("Dime la probabilidad del estado de tu interés "))
  probs<-as.numeric(probs)
  veces<-readline(prompt=paste("Dime el número de veces que quieres simular el escenario "))
  veces<-as.numeric(veces)
  (factorial(veces)/(factorial(x)*factorial(veces-x)))*(probs^x)*((1-probs)^y)
}
binomial.prob(8,0)
binomial.prob(10,0)
binomial.prob(7,3)
binomial.prob(3,7)


#----galaxy_connector----
install.packages("GalaxyConnector")
sxsx
#----single-cell----
#install packages 
BiocManager::install('SingleCellExperiment')
BiocManager::install(c('scater', 'scran', 'uwot'))
BiocManager::install('AnnotationHub')
BiocManager::install('scRNAseq')
BiocManager::install('BiocFileCache')
corral#tal vez sea de cran

library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(AnnotationHub)
library(scRNAseq)
library(BiocFileCache)

BiocManager::install(c("DADA2","microbiome","phyloseq"))
library(dada2)
library(microbiome)
library(phyloseq)

#----diversities----
shan<-function(vc){
  ff<-c()
  ll<-c()
  for (i in 1:length(vc)) {
    ff[i]<- vc[i]/sum(vc)
    ll[i]<- -(ff[i]*log(ff[i],2))
  }
  kk<-sum(ll)
  nom<- kk/log(length(vc),2)
  print(paste("la diversidad es", kk, "la diversidad normalizada es", nom))
}

species.vector<-c(10,6,4)
shan(species.vector)
spec.1<-rep(4,100)
spec.2<-c(20,1,0,0,0,0,0,0,0)
spec.3<-rnorm(100,10,1)
shan(spec.1)
shan(spec.2)
shan(spec.3)

#--- 
simpson.fun<-function(v){
  rr<-c()
  for (i in 1:length(v)) {
    rr[i]<-(v[i]/sum(v))^2
  }
  ll<-sum(rr)
  return(1-ll)
}
ss<-c(10,6)
simpson.fun(spec.3)
#---
simpson.inv<-function(v){
  return(1/simpson.fun(v))
}
simpson.inv(ss)

# thing for the final project, it is to assign taxa with a different data set for training, therefore had to download some othe trhings


"~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/"
"~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/"

path.1<-"~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/"
dada2:::makeTaxonomyFasta_RDP(file.path(ruta.2, "fungiLSU_train_012014_lsu_fixed_v2.fa"), file.path(ruta.2, "fungiLSU_taxid_012014.txt"),"~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDP_LSU_fixed_train_set_v2.fa",compress=FALSE)
dada2:::makeSpeciesFasta_RDP("~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/current_Fungi_unaligned.fa", "~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/rdp_species_assignment_LSU_v2.fa", compress=FALSE)





dada2:::makeTaxonomyFasta_RDP("~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/fungiLSU_train_012014_lsu_fixed_v2.fa", file.path(path.1, "fungiLSU_taxid_012014.txt"),"~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDP_LSU_fixed_train_set_v2.fa.gz",compress=TRUE)
dada2:::makeTaxonomyFasta_RDP("~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/fungiLSU_train_012014_lsu_fixed_v2.fa", file.path(path.1, "fungiLSU_taxid_012014.txt"),"~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDP_LSU_fixed_train_set_v2.fa.zip",compress=TRUE)
dada2:::makeSpeciesFasta_RDP("~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/current_Fungi_unaligned.fa", "~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/rdp_species_assignment_LSU_v2.fa.gz", compress=TRUE)
dada2:::makeSpeciesFasta_RDP("~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/current_Fungi_unaligned.fa", "~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/rdp_species_assignment_LSU_v2.fa.zip", compress=TRUE)


dada2:::makeTaxonomyFasta_RDP("/media/lauren/96BA-19E6/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/fungiLSU_train_012014_lsu_fixed_v2.fa", file.path(path, "fungiLSU_taxid_012014.txt"),"/media/lauren/96BA-19E6/Upload/RDP_LSU_fixed_train_set_v2.fa.gz",compress=TRUE)
dada2:::makeTaxonomyFasta_RDP("/media/lauren/96BA-19E6/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/fungiLSU_train_012014_lsu_fixed_v2.fa", file.path(path, "fungiLSU_taxid_012014.txt"),"/media/lauren/96BA-19E6/Upload/RDP_LSU_fixed_train_set_v2.fa.zip",compress=TRUE)
dada2:::makeSpeciesFasta_RDP("~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/current_Fungi_unaligned.fa", "~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/rdp_species_assignment_LSU_v2.fa.gz", compress=TRUE)
dada2:::makeSpeciesFasta_RDP("~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/current_Fungi_unaligned.fa", "~/R/Genomica_funcional/pro_fin_genfun/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata/rdp_species_assignment_LSU_v2.fa.zip", compress=TRUE)
  
library(ggplotly)
library(RCy3)
