#----Primer ejercicio----
library(msa)
library(Biostrings)
library(seqinr)
library(igraph)
##caragar archivo
insu<-readAAStringSet("Insulinas.fasta")
##alineamientos
al_1<-msa(insu, method = "ClustalW")
al_2<-msa(insu, method = "ClustalOmega")
##convertir en matriz de distancias
al_1<-msaConvert(al_1, type="seqinr::alignment")
al_2<-msaConvert(al_2, type="seqinr::alignment")
dist_1<-dist.alignment(al_1)
dist_2<-dist.alignment(al_2)
mat_1<-as.matrix(dist_1)
mat_2<-as.matrix(dist_2)
##hacerlo en red
red_1<-graph_from_adjacency_matrix(mat_1, mode = "undirected", weighted = T)
plot(red_1,edge.width=E(red_1)$weight*5, edge.color="blue")
#clusterizar
bb<-label.propagation.community(red_1, weights = E(red_1)$weight)
cc<-cluster_optimal(red_1, weights = E(red_1)$weight)
dd<-cluster_spinglass(red_1, weights = E(red_1)$weight)
#discusión
#en todos los métodos de clusterización me salió un solo grupo, esto se fue debido probablemente a que no puse un umbral 
#en el cual se denomine como que hay una conexión o no, dado que el peso de todas termina atrayéndolas al cluster, parece ser que es suficiente el 
#peso como para generar continuamente este mismo cluste una y otra vez; aunque se puede ver que hay un nodo en particular que se aleja de los demás 
#cuando se grafica, quizás algún método más riguroso pudiera generar otro cluster

#----Segundo ejercicio----
library(igraphdata)
#cargar datos
data("yeast")
#verificar unas cosas
class(yeast)
is.directed(yeast)
is_weighted(yeast)
#poner los nombres de los nodos más conectados
names(sort(degree(yeast), decreasing = T)[1:10])
#gráfica de distribución
hist(degree(yeast), col = c(2:15))
#diámetro de red
diameter(yeast)
#coeficiente de clusterización de los 20 más conectados
aa<-sort(degree(yeast), decreasing = T)[1:20]
transitivity(yeast, "local",aa)
#Encontrar los que tienen coeficiente de clusterización de 1
ll<-which(transitivity(yeast,"local",isolates = c("zero")) == 1) 
V(yeast)$name[ll]
#el hecho que tengan coeficiente de 1 significa que es una red muy conectada, que tiene redundancia y que probablemete es de
#mundo pequeño o ultrapequeño. Estos nodos están conectado de modo que entre ellos todas susconexiones conectan con otro del grupo
#El porcentaje de conexiones respecto al total.**** EN ESTA ASUMO QUE TE REFIERES A LOS 20 MÁS CONECTADOS, NO ESTABAS EN ZOOM, ENTONCE LO ASUMÍ******
sum(aa)/length(E(yeast))
#promedio de conexiones
mean(degree(yeast))
#nodos importantes con 3 métodos distintos
sort(degree(yeast), decreasing = T)[1:10]
sort(eccentricity(yeast), decreasing = F)[1:10]
sort(closeness(yeast), decreasing = F)[1:10]
#más alejadas

#eliminar y medir
diam<-function(red){
  for (i in 1:100) {
    red<-delete.vertices(red,V(red)[sample(length(V(red)), 1)])
    print(mean())
  }
}
diam(yeast)


#----Tercer ejercicio----
#cargar tabla
ff<-read.csv("Red de co-gustos - Hoja 1-3.csv")
#Hacer la red
rownames(ff)<-ff[,1]
ff<-ff[,-1]
ff<-t(as.matrix(ff))
nn<-cor(ff, method = "pearson")
nn<-(nn+1)/2
diag(nn)<-rep(0,15)
rr<-graph_from_adjacency_matrix(nn, mode="undirected", weighted=T)
#aquí estoy consciente de que no puse un umbral para separar la red, solo quiero ver qué onda
#quiénes están juntos
bb<-label.propagation.community(rr, weights = E(rr)$weight)
dd<-cluster_spinglass(rr, weights = E(rr)$weight)#si no hacemos lo del umbral. este método de cluster nos da dos grupos de alumnos, solo se 
#debe de correr el objeto y listo, se ve quiénes son
#distanacias y heatmap
distances(rr)
heatmap(distances(rr))
#graficar con tres layouts
plot(rr, layout=layout.circle, main="círculo")
plot(rr, layout=layout.random, main="random")
plot(rr, layout=layout.fruchterman.reingold, main="fruchterman.reingold")
#graficar con medidas de centralidad
plot(rr, vertex.size=closeness(rr, mode = "all")*100,
     edge.width=E(rr)$weight*3, edge.color= "blue")#con closeness
plot(rr, vertex.size=degree(rr, mode = "all")+4)#con degree, sé que van a ser todos iguales porque tienen mismo # de conexiones
plot(rr, vertex.size=betweenness(rr)+7)#con betweeness
