#Ejercicios de genómica

#----prácticas_de_igraph----
library(igraph)
#1.- Calculate the degree,average degree, plot the degree distribution, obtain the adjacency matrix from the
#next networks
g1<-barabasi.game(100,directed = FALSE)
g2<-random.graph.game(100,0.5)
g3<-sample_smallworld(1,100,p=0.2,nei=3)
#degree
degree(g1)
degree(g2)
degree(g3)
#average degree
mean(degree(g1))
mean(degree(g2))
mean(degree(g3))
#plot degree distribution
plot(degree.distribution(g1))
plot(degree.distribution(g2))
plot(degree.distribution(g3))
#adjancy matrices
as.matrix(get.adjacency(g1))
as.matrix(get.adjacency(g2))
as.matrix(get.adjacency(g3))

#2.- Identify the vertex whose degree is the maximum
sort(degree(g1), decreasing = T)[1]
sort(degree(g2), decreasing = T)[1]
sort(degree(g3), decreasing = T)[1]

#3.- Find the 10 most conneceted vertices
sort(degree(g1), decreasing = T)[1:10]
sort(degree(g2), decreasing = T)[1:10]
sort(degree(g3), decreasing = T)[1:10]

#4.- Build a network from vertex = people form this classroom, rule of connection= two people are connected
#if they were born in a neighbor state or in the same state. Use colors, names and sizes that reflects the
#totpological properties of the network.
#cannot make, but would be interesting


#----prácticas_de_distances_clustering_coefficient----
library(igraph)
#1.- Using the networks from slide number calculate the distances, mean of distances,the path from the
#most distant nodes and the clustering coefficient
#do not know which slides it is referring, but still, it would be like this
#distance
distances()
#mean distance
mean_distance()
#path from most distant nodes
get_diameter(red_1)#shows the path directly
#clustering coefficient
transitivity()#probability that the adjacent vertices of a vertex are connected

#2.- Using the next networks calculate the distances, mean of distances, the path from the most distant
#nodes and the clustering coefficient
g1<-barabasi.game(100,directed = FALSE)
g2<-random.graph.game(100,0.20)
g3<-sample_smallworld(1,100,p=0.2,nei=3)
#distances
distances(g1)
distances(g2)
distances(g3)
#mean distances
mean_distance(g1)
mean_distance(g2)
mean_distance(g3)
#paths from most distant nodes 
get_diameter(g1)
get_diameter(g2)
get_diameter(g3)
#clustering coefficient
transitivity(g1)
transitivity(g2)
transitivity(g3)

#3.-5.Build a network from vertex = people form this room, rule of connection= two people are connected if they
#were born in the same state. calculate the distances, mean of distances, the path from the most distant nodes
#and the clustering coefficient
#cannot do, but would be nice 


#----prácticas_de_free_scale_property*incomplete*----
library(igraph)
#1.- From the next networks, plot the degree distribution
g10<-barabasi.game(10,directed = FALSE)
g100<-barabasi.game(100,directed = FALSE)
g1K<-barabasi.game(1000,directed = FALSE)
g2<-random.graph.game(1000,0.20)
g3<-sample_smallworld(1,1000,p=0.2,nei=3)
g4<-make_graph("Zachary")
#degree distribution
plot(degree_distribution(g10))
plot(degree_distribution(g100))
plot(degree_distribution(g1K))
plot(degree_distribution(g2))
plot(degree_distribution(g3))
plot(degree_distribution(g4))

#2.- Find the median, mean and boxplots of these distributions.
#median
median(degree_distribution(g10))
median(degree_distribution(g100))
median(degree_distribution(g1K))
median(degree_distribution(g2))
median(degree_distribution(g3))
median(degree_distribution(g4))
#mean
mean(degree_distribution(g10))
mean(degree_distribution(g100))
mean(degree_distribution(g1K))
mean(degree_distribution(g2))
mean(degree_distribution(g3))
mean(degree_distribution(g4))
#boxplots
boxplot(degree_distribution(g10))
boxplot(degree_distribution(g100))
boxplot(degree_distribution(g1K))
boxplot(degree_distribution(g2))
boxplot(degree_distribution(g3))
boxplot(degree_distribution(g4))

#3.- Fit a power law to these distributions. Discuss your results
fit_power_law(degree(g10)+1, 10)
fit_power_law(degree(g100)+1, 10)
fit_power_law(degree(g1K)+1, 10)
fit_power_law(degree(g2)+1, 10)
fit_power_law(degree(g3)+1, 10)
fit_power_law(degree(g4)+1, 10)
#being honest, i don't have a single clue what this means

#4.- Fit a power law distribution to the classroom network

#I don't really know whta this means, must read before making anything


#----prátcias_network_robustness----
library(igraph)
#1.- Using the next networks, remove randomly 1/100 nodes. Calculate the mean of the distance and the
#diameter. Plot the network. Repeat this process 10 times for each network.
g100<-barabasi.game(100,directed = FALSE)
g1K<-barabasi.game(1000,directed = FALSE)
g2<-random.graph.game(1000,0.20)
g3<-sample_smallworld(1,1000,p=0.2,nei=3)
#cannot be bothered doing it manually, made a function, decided to omit the plot part 'casue it would 
#take forever in my computer
rem.rd.1<-function(z){
  for (i in 1:10) {
    z<-delete.vertices(z,sample(1:length(V(z)),1))
    print(paste("the mean distance for",i, "deletions is", mean_distance(z)))
    print(paste("the diameter for",i, "deletions is", diameter(z)))
    #print(plot(z)) DO NOT RUN THIS LINE UNLESS YOU HAVE A *VERY* GOOD COMPUTER, UNLESS YOU WANT 
    #TO WASTE YOUR TIME
  }
}

rem.rd.1(g100)
rem.rd.1(g1K)
rem.rd.1(g2)
rem.rd.1(g3)
#2.- Using the same networks, remove the 10 most connected nodes. Calculate the mean of the distance and
#the diameter. Plot the network. Repeat this process 10 times for each network.
#again, couldn't be bothered doing it manually, just modified previous function 
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


#----librería_red_de_amigos----
library(igraph)#cargar la lirería
fr<-read.csv("RedDeAmistadesGF2021Nombres - Hoja 1.csv")#así se llama el archivo que descargué
rownames(fr)<-fr[,1]
fr_mat<-as.matrix(fr[,-1])
diag(fr_mat)<-rep(0,19)#eliminar los datos que estaban en la tabla que nos harían ruido
red_1<-graph_from_adjacency_matrix(fr_mat, mode= "directed")#hacer una red pesada co la tabla ya modificada
red_1


#----Exercise figures first time----
#networks
library(igraph)
estrella<-make_star(10, mode = "undirected")
anillo<-make_ring(10, directed = F)
above_right<-make_empty_graph(n=10, directed=F)
above_right<-add.edges(above_right, c(1,10, 1,9, 1,2, 1,8, 1,4, 2,10, 2,9, 
                                      2,7, 2,6, 9,10, 9,7, 9,8, 9,4, 7,4,
                                      7,5, 7,3, 3,5, 3,8, 6,4, 6,8, 6,5, 7,10,
                                      8,5, 8,4, 5,4))
down_right<-make_empty_graph(n=10, directed=F)
down_right<-add.edges(down_right, c(4,3, 3,1, 8,1, 9,2, 2,1, 6,1, 5,1,
                                    7,5, 10,5))
#1.-Star
##Degree
degree(estrella)
#nodo 1---- 9
#resto----- 1
##Av. degree
mean(degree(estrella))
#1.8
##Degree distribution
plot(degree.distribution(estrella))
#0.9---1 connection
#0.1---9 connections

#2.-Ring
#Degree
degree(anillo)
#todos 2
#Av. degree
mean(degree(anillo))
#2
#Degree distribution
plot(degree.distribution(anillo))
#100%---2

#3.-Above Right
#Degree
degree(above_right)
#1,2,5----5
#3------3
#4,7,8,9------6
#6,10------4
#av. degree
mean(degree(above_right))
#5
#degree distribution
plot(degree.distribution(above_right))
#0.3----5
#0.1----3
#0.4----6
#0.2----4

#4.-down right
#Degree
degree(down_right)
#1------5
#2,3------2
#4,6,7,8,9,10------1
#5------3
#Av. degree
mean(degree(down_right))
#1.8
#Degree distribution
plot(degree.distribution(down_right))
#0.1----5
#0.2----2
#0.6----1
#0.1----3


#----igrpah_first_exercises----
library(igraph)
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

#----Exercises_with_classmates----
library(igraph)#cargar la lirería
fr<-read.csv("RedDeAmistadesGF2021Nombres - Hoja 1.csv")#así se llama el archivo que descargué
rownames(fr)<-fr[,1]
fr_mat<-as.matrix(fr[,-1])
diag(fr_mat)<-rep(0,19)#eliminar los datos que estaban en la tabla que nos harían ruido
red_1<-graph_from_adjacency_matrix(fr_mat, mode= "directed")#hacer una red pesada co la tabla ya modificada
red_1
###
plot(red_1)
degree(red_1, mode="in")
plot(cluster_spinglass(red_1),red_1)#plotear con un método de clusterización, en este caso es el de cluster_spinglass
hist(degree(red_1,mode="in"), col=c(1:length(hist(degree(red_1, mode="in")))))
hist(degree(red_1,mode="out"), col=c(1:length(hist(degree(red_1, mode="out")))))

#----Second_ex_friend_network----
library(igraph)
#the network is tyhe one from the previous work

#Distances

shortest.paths(red_1)#shortest path between a node and the others
distances(red_1, v= "ANA", to="ENOE")#may select a single node to other
diameter(red_1)#diametro de la red, checar lo de los argumentos que se tienen de modo estándar
mean_distance(red_1)#mean distance of the distances
plot(red_1,layout=layout.fruchterman.reingold(red_1), #plotear la red con el layout que se quiere
     vertex.size=5, edge.arrow.size=.2)
###
g <- barabasi.game(1000, power=1)#hacer una red basada en el modelo de barabasi estocástico
plot(g,layout=layout.fruchterman.reingold(g), 
     vertex.size=5, edge.arrow.size=.2)
###
mat.1<-shortest.paths(red_1)#apparently seems to be a matrix by default, matrix of all the shortest paths between all the nodes
heatmap(mat.1)
shortest_paths(red_1,1,6)
all_shortest_paths(red_1, 1,2:6)




#----exercise_find_solution-----
#generar función rudimentaria para encontrar valor más grnade de un vector
a<-c(13,21,54,42,3,54)
a<-c(13,21,54,42,3,54,100)
for (i in 1:length(a)) {
  if(sum(a[i]>=a[1:length(a)])==(length(a))){
    print(paste(a[i], "es el más grande"))
  }
}

#----ejercicios_robustness_y_otras_cosas-----
alt_h<-c(1.70, 1.73, 1.82, 1.76, 1.75, 1.76)#alturas hombres
alt_m<-c(1.63, 1.59, 1.65, 1.69, 1.70, 1.60, 1.68, 1.64, 1.59, 1.58, 1.63)#alturas mujeres
hist(alt_h)#representaciones en histograma
hist(alt_m)
sd(alt_h)#desviaciones estándar
sd(alt_m)
mean(alt_h)#medias
mean(alt_m)

red_1#cargar otra vez la red de amigos del salón

barplot(sort(degree(red_1, mode = "in")), beside = T, col = rainbow(19))#plotear con colores basados en el degree de entrada
barplot(sort(degree(red_1, mode = "out")), beside = T, col = rainbow(19))#lo mismo que arriba pero de salida

plot(sort(degree(red_1, mode = "in")))#mismo que las de arriba pero con los puntos en el gráfico en lugar de barras
plot(sort(degree(red_1, mode = "out")))

sort(degree(red_1, mode = "in"), decreasing = T)[1:3]#los de mayor degree de entrada
sort(degree(red_1, mode = "out"), decreasing = T)[1:3]#los de mayor degree de salida
####

red_1["JORGE"]#obtener en qué nodos conecta el nodo deseado con otros (el 1) 
as.vector(red_1["JORGE"])#mismo que línea de arriba, pero lo paso a un vector
aa<-as.vector(names(which(red_1["JORGE"]==1)))#obtener nombre de los nodos que el nodo de interés conecta
ma.1<-as.matrix(red_1[aa])#obtener las conexiones que poseen los nodos que fueron conectados a partir del de interés
mean(rowSums(ma.1))#promedio de conexiones que poseen los nodos a los que conecta el nodo de interés
#

#------robustness

g100<-barabasi.game(100,directed = FALSE)
g1K<-barabasi.game(1000,directed = FALSE)
g2<-random.graph.game(1000,0.20)
g3<-sample_smallworld(1,1000,p=0.2,nei=3)
##
#función para remover un nodo a la vez de una red deseada y calcular l
#a distancia promedio de la red y su diámetro
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
#función para remover un nodo a la vez de una red deseada y calcular l
#a distancia promedio de la red y su diámetro
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


#----medidas_de_centralidad----

library(igraph)

a<-barabasi.game(100, directed = F)
degree(a)
sort(degree(a), decreasing = T)
eccentricity(a)#medida de excentricidad
sort(eccentricity(a),decreasing = F)
closeness(a)#medida de cercanía
betweenness(a)#medida de betweenness


sort(degree(a), decreasing = T)
sort(eccentricity(a),decreasing = F)
sort(closeness(a),decreasing = T)
sort(betweenness(a),decreasing = T)
eigen_centrality(a)#centralidad de eigenvalores
page.rank(a)
#with my network
#con la red que se ha usado desde casi el primer ejercicio
sort(degree(red_1, mode = "in"), decreasing = T)
sort(degree(red_1, mode = "out"), decreasing = T)
sort(eccentricity(red_1, mode = "in"),decreasing = F)
sort(eccentricity(red_1, mode = "out"),decreasing = F)
sort(closeness(red_1, mode = "in"),decreasing = T)
sort(closeness(red_1, mode = "out"),decreasing = T)
sort(betweenness(red_1, directed = T),decreasing = T)
eigen_centrality(red_1, directed = T)#creo que es en relación de los nodos con más conexiones, esos nodos qué tan conectados con los demás
page.rank(red_1)

#----generar_red_de_matriz_adyacencia----
library(igraph)
mm<-matrix(c(0,1,1,1,0,0,1,0,0), ncol = 3, nrow = 3)#generar una matriz de adyacencia
rr<-graph_from_adjacency_matrix(mm, mode = "undirected")#generar la red a partir de la matriz de adyacencia
V(rr)$name<-letters[1:3]#ponerle nombre a los modos
eigen_centrality(rr)[1]#centralidad de eigen
eigen(mm)#valores de eigen 
#----libreria_igraphdata----
library(igraphdata)
data(package="igraphdata")
data(yeast)
yeast
#----plotear_clusters----
#otra vez la misma red de siempre para cargar
#centrado en plotear las agrupaciones
red_1
aa<-edge.betweenness.community(red_1,directed = T, modularity = NULL)
bb<-label.propagation.community(red_1, weights = NULL)
cc<-cluster_optimal(red_1, weights = NULL)
dd<-cluster_spinglass(red_1)
plot(aa, red_1)
plot(bb, red_1)
plot(cc, red_1)
plot(dd, red_1)
plot(leading.eigenvector.community(red_1),red_1)
#----idk----
#try the values after applying log to large numbers
log(10000100000)
log(1000000)
log(log(10000100000))
log(log(1000000))
#make random networks
aa<-barabasi.game(100, directed = F)
bb<-barabasi.game(1000, directed = F)
cc<-random.graph.game(100, 1/10)
dd<-random.graph.game(1000, 1/10)
ee<-sample_smallworld(1, 100, 5, 0.1, 0.05)
#basic info from networks
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
red_1#hasta aquí es lo mismo de la red, pero modificada un poco
get.adjlist(red_1, mode = "out")#sacar la lista de adyacencias, conexiones al menos unilaterales
ff<-as.matrix(as_edgelist(red_1, names = T))#lista con los nodos de salida y los receptores
write.csv(ff, "prueba.csv")#generar el archivo del producto de la línea de arriba

#----red_de_gustos----
#redes de gustos 
gg<-read.csv("Red de co-gustos - Hoja 1.csv")#cargar la red
rownames(gg)<-gg[,1]#poner los nombres a cada renglón
gg<-gg[,-1]
gg<-t(as.matrix(gg))#transponer
#sacar diferentes valores de correlaciones
cor(gg, method = "pearson")
cor(gg, method = "kendall")
cor(gg, method = "spearman")
library(pheatmap)
pheatmap(cor(gg, method = "pearson"))#generar el heatmap
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
plot(cluster_spinglass(rr, weights = E(rr)$weight),rr)
##Formar red pesada
library(WGCNA)
library(igraphinshiny)
library(igraph)
rr<-graph_from_adjacency_matrix(mm, mode="undirected", weighted=T)
rr
plot(rr,edge.width=E(rr)$weight*2.5, edge.color="blue")
write.csv(mm, "pesada_1.csv")
hh<-matrix(as_edgelist(rr, names = T), ncol = 2)


#----boolean_networks_examples----
library(BoolNet)
red_3<-loadNetwork("Primera_red.txt")#cargar una red booleana 
plotNetworkWiring(red_3)#plotear la red
atrac<-getAttractors(red_3)#conseguir los atractores
plotAttractors(atrac)#plotear los atractores
plotSequence(red_3, includeAttractorStates = "all", startState = rep(1,2))#not working
getPathToAttractor(red_3, state = c(0,1,1), includeAttractorStates = "all")

red_pina<-loadNetwork("new_2.txt")
plotNetworkWiring(red_pina)
atrac<-getAttractors(red_3)
plotAttractors(atrac)

#----tutorial_WGCNA*incompleto*----
library(igraph)
library(WGCNA)
#cargar los datos que ya vienen en el tutorial
clinical_traits<-read.csv("FemaleLiver-Data/ClinicalTraits.csv")
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

#----red_booleana_cellcycle----

library(BoolNet)
data(cellcycle)
plotStateGraph(getAttractors(cellcycle))

#----Ejercicio_red_booleana----
library(igraph)
library(BoolNet)
lac<-loadNetwork("lac_bool.txt")#cargar la red
lac
plotNetworkWiring(lac)#plotear la red
ff<-generateState(lac, c("M"=1, "P"=1))#generar un estado artificial
getAttractors(lac, startStates = list(ff))#generar los atractores a partir del estado artificial generado

#----ejercicios de langebio----
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

            
#----librerias_uso_de_datos----
BiocManager::install("rhdf5")
devtools::install_github("pachterlab/sleuth")
library(devtools)
library(rhdf5)
library(sleuth)
library(edgeR)
library(ensembldb)

setwd("C:/Users/fotgo/OneDrive/Documentos/R/Genomica_funcional/")

#----Calcular_prueba_estadística_binomial----

#10 soles, 5 y 5, 3 y 7, 1 y 9

binomial.prob<-function(x,y){
  probs<-readline(prompt=paste("Dime la probabilidad del estado de tu interés "))
  probs<-as.numeric(probs)
  veces<-readline(prompt=paste("Dime el número de veces que quieres simular el escenario "))
  veces<-as.numeric(veces)
  (factorial(veces)/(factorial(x)*factorial(veces-x)))*(probs^x)*((1-probs)^y)
}
binomial.prob(8,0)
ll(10,0)
ll(7,3)
ll(3,7)

