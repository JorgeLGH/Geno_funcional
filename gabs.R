#### EXAMEN PRIMER PARCIAL #####

#### PROGRAMA 1 ####
#para esta parte necesitamos el paquete: igraph y igraphdata
data("yeast")

#1. degree
class (yeast) #efectivamente es un objeto tipo igraph 

degree (yeast) -> dy

sort(dy, decreasing = TRUE) #Odenamos de mayor a menor los degree

sort(dy, decreasing = TRUE) [1:10] #nos indica cuales son los 10 mas
#los nombres de cada proteína vienen junto con el degree

####################################################

#2. La gráfica de la distribución de conectividades.
hist (dy,col = "pink", main = paste("Distribución de conectividades de yeast"), xlab = "degree")
#como se observa en la grafica de distribución, la mayoria de los nodos tienen poquitas conexiones
#y hay muy pocos nodos con muchisimas conexiones
plot_dd <- function(n){
  plot (degree.distribution(n))
}

plot_dd(yeast)
#Al igual podriamos ver la distribucion con esta grafica y es lo mismo
#muchos nodos con pocas conexiones y muy pocos nodos con muchas conexiones

#pareciera ser una red free-scale esta red de proteínas de alguna manera 
#esta protegida de ataques dirigidos, ya que es más probable que el 
#ataque caiga en un nodo con muy poquitas conexiones y por lo cual no afecte a la red

##########################################################

#3. El diámetro de la red
diameter(yeast)
#El resultado es 15, lo que indica que es la maxima 
#trayectoria que existe entre dos nodos de la red

###########################################################

#4.  coeficiente de clusterización cada una 
#de las 20 proteínas más conectadas.

sort(dy, decreasing = TRUE) [1:20] -> y20 
#estas son las 20 proteinas mas conectadas,lo que necesitamos saber
#es cual es su numero en la red para saber que numero de la red es 
#para saber que numero del coeficiente es

for (i in y20)
{
  which (dy == i) -> y
  print (transitivity(yeast, type = "local") [y])
}

# para saber cual es el coef. de las 20 mas conectadas lo que hice fue 
#un ciclo for donde va a tomar cada uno de los degree de las 20 más conectadas
#y lo va comprar con el degree original, nos va a sacar una lista donde van a 
#venir la posicion de cada nodo y la vaa buscar con la funcion de transitivity
#por lo cual nos dara el coeficiente de todos los 20 nodos más conectados

###################################################

#5. proteínas con coeficiente de clusterización de 1
which(transitivity(yeast,type= c("local"), vids = NULL,
                   isolates = c("zero")) == 1) 
#nos esta arrogando todos los nodos que tienen un coeficiente igual a 1
#son varios y lo que significa es que todos esos nodos estan interaccionando entre si
#y es una red muy conectada, puede ser que sea una red small world

######################################################

#6. porcentaje de conexiones respecto al total
length(yeast)
sum(dy) #nos indica cual es el todal de todas las conexiones
#las cuales son el 23710 y son 10 nodos



#####################################################

#7. El promedio de conectividades
mean (dy) -> promedio
promedio #el promedio es de 9, eso quiere decir que en promedio
#los nodos de esta red de proteínas tienen 9 conexiones

#####################################################

#8. Encuentre los caminos entre una pareja de las 10 proteínas más conectadas
sort(dy, decreasing = TRUE) [1:10]

#Voy a seleccionar a la proteina YPR110C (con degree de 118) y a la YPL131W (degree de 115)
which(dy == 118) #su nodo es el 286
which(dy == 115) #su nodo es el 698

distances(yeast, v = 118, to = 698, mode = "all") #al parecer solo hay camino entre estas dos proteinas

########################################################

#9. Seleccione 1 proteína al azar y la la elimine de la red, que  calcule el promedio de las distancias después de quitar la proteína. 
#Hacer esto 10 veces más y compare los promedios.

diame <- function(red) {
  for (i in 1:10) 
  {
    rf <- delete.vertices(red, sample (1:(10-i),1))
    print (mean_distance(rf))
    diam <- diameter(rf)
    print(diam)
  }
}

diame (yeast)
#lo que nos imprime primero es el diametro despues de quitar un nodo
#abajo esta el diamtero, y lo esta comparando con el priomedio

plot(1:10,diame(yeast), type= "l") #y lo que observamos en la grafica
#es que ladistancia va a aumentar cuando quitemos nodos 

##################################################
#10. Discusión de resultados

#lo que obtenemos en esta red, a lo largo de los ejercicios primero es que
#es una red que tiene muchos nodos con pocas conexiones y muy pocos nodos con muchas
#conexiones, caracteristica de redes free-scale, tambien tiene muchos nodos con 
#coef. de cluster de igual a 1, eso quiere decir que sus nodos estan muy conectados entre si
#y es una caracteristica de redes small world, al quitarle nodos se observa lo esperado
#que la distancia aumenta ya que esos nodos deben buscar otros caminos para seguir conectados.
#Es importante recordar que esta es una red biologica de proteinas, como sabemos las proteinas son muy
#importantes para la supervivencia de un organismo y que no se pueden dar el lujo de perder por lo cual
#deben estar protegidas, de manera que esta diseñada para el supuesto de un ataque dirigido no afecte a las proteinas
#más conectadas y no afecte su integridad o en el caso de ser atacadas, puedan encontrar caminos alternos para seguir conectados
#y funcionando. 


#### PROGRAMA 2 ####
#para esta parte necesitamos el paquete: msa y seqinr

readAAStringSet("K:/GENOMICA_FUNCIONAL_GABRIELA/Insulinas.fasta") -> insulinas
insulinas

#Nota: al aparecer mi version de r no tiene la funcion msaConvert
#a pesar de tener las librerias msa y seqinr, por lo cual escribire 
#como se hace,espero que corra :(

a_insulinas <- msa(insulinas) #alineamiento
a_convertido <- msaConvert(a_insulinas, type = "seqinr::alignment")
a_distancia <- dist.alignment (a_convertido, "identity") #ya transformado, lo ponemos en un alineamiento de distancias
a_matriz <- as.matrix(a_distancia) #lo pasamos a matriz, para que se puede leer despues en la red

red_insulinas <- graph_from_adjacency_matrix(a_matriz) #para pasarla a una red

#por ultimo, ya que tenemos la red, nos piden hacer un cluster 
#esto se haria de la siguiente manera
cluster_A <- cluster_edge_betweenness(red_insulinas)
plot(cluster_A, red_insulinas) #nos visualiza los grupos que secrean
# en cuestion de como los hizo este metodo de clusterizacion 
