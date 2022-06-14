#Tarea 01

#1.-A partir de las redes de la figura 1. Calcula con igraph, las siguientes propiedades:

library(igraph)
library(bazar)#necesaria para la función

estrella<-make_star(10, mode = "undirected")
anillo<-make_ring(10, directed = F)
up_right<-make_empty_graph(n=10, directed=F)
up_right<-add.edges(up_right, c(1,2, 1,4, 1,9, 1,10, 1,8, 10,2, 10,9, 10,7, 2,9, 2,7, 2,6,
                                9,7, 9,8, 9,4, 7,3, 7,4, 7,5, 3,8, 3,5, 8,4, 8,5, 8,6, 
                                5,6, 5,4, 4,6))
down_right<-make_empty_graph(n=10, directed = F)
down_right<-add.edges(down_right, c(1,6, 1,3, 1,8, 1,2, 1,5, 3,4, 2,9, 5,7, 5,10))

basic_info<-function(network){
  print(paste("la distribución de conectaivides es de", 
              concat(degree.distribution(network), sep = ", ")))
  print(paste("el nodo más conectado es/son ", concat(V(network)[which(degree(network)==max(degree(network)))]) ,
              "y tiene/n ", concat(sort(rowSums(as.matrix(get.adjacency(network))), 
                                                decreasing = T)[1]), "conexiones"))
  print(paste("el diámetro de la red es ", diameter(network)))
  print(distances(network))
  print(heatmap(distances(network)))
  for (i in 1:length(V(network))) {
    if(sum(length(neighbors(network, i)))==1){
      print(paste("el vecino de ", i, "es", neighbors(network, i)))
    }else print(paste("los vecinos de ", i, "son", concat(neighbors(network, i), sep = ", ")))
  }
}

basic_info(down_right)
basic_info(estrella)
basic_info(anillo)
basic_info(up_right)

#2.- Elabora un programa en R que utilice un ciclo for para a partir del vector v siguiente

v<-sample(100)

for (i in 1:length(v)){
  if((i%%2)==1){
    print(i^2)
  }
}

#3.- Elabora un programa en R que a partir del archivo de amistades del grupo.

#a)
fr<-read.csv("RedDeAmistadesGF2021Nombres - Hoja 1.csv")
rownames(fr)<-fr[,1]
fr_mat<-as.matrix(fr[,-1])
diag(fr_mat)<-rep(0,19)
red_1<-graph_from_adjacency_matrix(fr_mat, mode= "directed")

#b)
red_1["JORGE"]
aa<-as.vector(names(which(red_1["JORGE"]==1)))
aa

#c)
bb<-as.vector(names(which(get.adjacency(red_1)[,8]==1)))
bb

#d)
for (i in 1:length(aa)) {
  print(paste("hola amig@", aa[i]))
}

#e)
mean(rowSums(as.matrix(get.adjacency(red_1))))

#f)
lista<-c(4,5,10,13,15,17,19)
g<-0
for (i in 1:length(aa)) {
  ss<-sum(red_1[aa[i]][lista])
  g[i]<-ss
  h<-sum(g)
}
h
h/((sum(red_1["JORGE"]))*(sum(red_1["JORGE"])-1))

#4.- Utiliza la red del club de Karate de Zachary (investiga cmo puedes generarla en igraph)

library(igraphdata)
data(package="igraphdata")
data(karate)
karate

#a)
sort(degree(karate), decreasing = T)[c(1:3)]
#b)
diameter(karate)
#c)
degree.distribution(karate)
#d)
get.adjacency(karate)
#e)
plot(karate, layout=layout_nicely, vertex.size=degree(karate, V(karate))*3+2,
     vertex.label.dist=0.5, edge.arrow.size=.5)
