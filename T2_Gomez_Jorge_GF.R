library(igraph)
library(BoolNet)
library(QCA)

#----Ejercicio 1----
sink("red_1")
cat("targets, factors\n")
cat("A, A & C | ! B \n")
cat("B, C | ! D\n")
cat("C, B | ! D\n")
cat("D, A | ! C\n")
sink() #con esto se genera la red

red_1 <- loadNetwork("red_1") #cargar la red que hice en un objeto

#tabla de verdad
#A(A&C,B)= (T,F) | ! T =F 
#A(A&C,B)= (T,T) | ! T =T
#A(A&C,B)= (F,T) | ! T =F 
#A(A&C,B)= (F,F) | ! T =F 
#A(A&C,B)= (T,F) | ! F =T 
#A(A&C,B)= (T,T) | ! F =T
#A(A&C,B)= (F,T) | ! F =T 
#A(A&C,B)= (F,F) | ! F =T 
#B(c,D)= T | ! F = T
#B(c,D)= F | ! F = T
#B(c,D)= T | ! T = T
#B(c,D)= F | ! T = F
#C(B,D)= T | ! F = T
#C(B,D)= F | ! F = T
#C(B,D)= T | ! T = T
#C(B,D)= F | ! T = F
#D(A,C)= T | ! F = T
#D(A,C)= F | ! F = T
#D(A,C)= T | ! T = T
#D(A,C)= F | ! T = F

getTransitionTable(getAttractors(red))#aquí están todos los posibles estados y su siguiente estado

atractores <- getAttractors(red) #aquí se obtienen atractores

plotAttractors(atractores) #plot de los atractores

#----Ejercicio 2----
sink("red_2")
cat("targets, factors\n")
cat("CIS, INS | CADPR & CGMP\n")
cat("CGMP, GC\n")
cat("CADPR, ADPRC\n")
cat("INS, PLC\n")
cat("GC, NO\n")
cat("ADPRC, NO\n")
cat("PLC, CA\n")
cat("NO, NOS\n")
cat("NOS, CA\n")
cat("CA, CALM | CIS | ! CATPASE\n")
cat("CATPASE, CA\n")
cat("CALM, ! DEPOLAR\n")
cat("KEV, CA\n")
cat("HATPASE, ! CA\n")
cat("ANION, CA\n")
cat("CLOSURE, ANION | KOUT & KAP\n")
cat("KOUT, DEPOLAR\n")
cat("KAP, DEPOLAR\n")
cat("DEPOLAR, KEV & ! HATPASE | ! KOUT | ANION\n")
sink()#generar la red, está todo en mayúsculas porque así no me confundo

red_2 <- loadNetwork("red_2")#cargar la red

truthTableToSymbolic(red_2)#aquí viene una representación de las reglas

atractores_2 <- getAttractors(red_2) #encontré un atractor con las reglas que establecí
#Es un atractor cíclico, de manera que probablemente el proceso al cual está asociado requiere de la continua activación
#e inhibición de algunos componentes para poder llevar a cabo su correcto funcionamiento. El sentido biológico es que ciertas enzimas van a estar
#permietiendo la expresión y funcionamiento de esta red, pero se autoregula de modo que cambia de estados de manera cíclica

#----Ejercicio 3----
sink("red_3")
cat("targets, factors\n")
cat("D, D\n")
cat("TWI, TWI | D | FOG\n")
cat("SNA, TWI | D\n")
cat("FOG, SNA | TWI\n")
cat("B, SNA\n")
cat("FOGR, FOG\n")
cat("CTA, FOGR\n")
cat("A, TWI | A\n")
cat("CSK, CTA\n")
cat("SRC, ! CSK\n")
cat("GAP, SRC\n")
cat("GEF, CTA | A\n")
cat("RHO, GEF & ! GAP & ! GDI\n")
cat("ROCK, RHO\n")
cat("MLCP, ! ROCK\n")
cat("MOE, MOE | ! MLCP | ROCK\n")
cat("GDI, ! MOE\n")
cat("ACTIN, B | MLC | MOE\n")
cat("MLC, ! MLCP | MLCK | ROCK\n")
cat("MLCK, MLCK\n")
sink()#generar red

red_3 <- loadNetwork("red_3")#cargar la red a un objeto

atractores_3 <- getAttractors(red_3)#los atractores

#pues dependiendo del estado inicial de las demás entradas el resultado que se dará
#las demás entradas limitan a qué atractor se acopla mejor el comportamiento inicial propuesto
#con la función de abajo podemos ver cómo se comporta 
stateTransition(red_3, c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1))


#----Ejercicio 4----
library(matlib)
mat_1 <- matrix(c(1,2,3,4),ncol = 2)#matriz con la que trabajaré

#no funciona, prometo que sé cómo hacerlo, pero pues no me da
eig <- function(matrix){
  mati <- matrix(c(0,0,0,0), ncol = 2)
  x<-1
  diag(mati)<-(x)
  ff<-matrix-mati
  rr<-ff[1,1]*ff[2,2]
  mm<-ff[1,2]*ff[2,1]
  jj<-rr-mm
  return(rr)
}

eig(mat_1)
eigen(mat_1)

#----Ejercicio 5----
#realizado con el de acceso GSE25724
#https://github.com/JorgeLGH/Ejercicio_5.git
#estoy tomando en cuenta que el script de GEO2R ya lo tienes, entonces solo procederé de los datos que ya obtuve para lo de las tablas en csv
deTable <- topTable(fit2, number=nrow(fit2), lfc=log2(2), p.value=0.1)
deTable
write.csv(deTable, "resultados.csv")
#----Ejercicio 6----
datos<-data.frame()
usuarios.1<-function(data_frame){
  pregunta<-readline(prompt=paste("¿Eres usuario nuevo? (Y/N) "))
  pregunta<-as.character(pregunta)
  if(pregunta=="y"|pregunta=="Y"|pregunta=="Yes"|pregunta=="YES"){
    nombre.n<-readline(prompt=paste("Designa tu nombre de usuario "))
    nombre.n<-as.character(nombre.n)
    contra.n<-readline(prompt=paste("Escribe tu contraseña "))
    contra.n<-as.character(contra.n)
    rep<-readline(prompt=paste("Repite tu contraseña "))
    rep<-as.character(rep)
    if(sum(rep==contra.n)==1){
      data_frame[nrow(data_frame)+1,1]<-nombre.n
      data_frame[nrow(data_frame),2]<-contra.n
      print(paste("Tu registro quedó completo "))
      return(data_frame)
    }else print(paste("No coincide la contraseña, inéntalo de nuevo"))
  }else
    respuesta<-readline(prompt=paste("Ingrese su nombre de usuario "))
  respuesta<-as.character(respuesta)
  if(sum(respuesta==data_frame[ ,1])==1){
    respuesta2<-readline(prompt=paste("Ingresa tu contraseña "))
    respuesta2<-as.character(respuesta2)
    if(sum(respuesta2==data_frame[ ,2])==1){
      respuesta3<-readline(prompt=paste("Confrima tu contraseña "))
      respuesta3<-as.character(respuesta3)
      if(respuesta3==respuesta2){
        print(paste("Bienvenido"))
      }else print(paste("No coincide la verificación de la contraseña"))
    }
  }else print(paste("Usuario incorrecto, intenténtalo de nuevo"))
}
datos<-usuarios.1(datos) #hay muchas cosas que mejorar, pero pues se hace mínimamente lo que se pide
#----Ejercicio 7----
#requiere otros objetos de sesiones pasadas, la verdad no lo hice por dejarlo al final
#y no pienso fingir que entiendo lo que se hace comentando exactamente lo que
#viene en el texto del tutorial
#----Ejercicio 8----
ejemplo<-data.frame()
ejemplo[c(1:4),1]<- c("Mario","Carmen", "Juan", "Daniela")
ejemplo[c(1:4),2]<- c(1,4,76,9)
ganadores<-function(data_frame){
  gana<-sample(ejemplo[,2], 1)
  print(paste("felicidades", ejemplo[which(ejemplo[,2]==gana),1], "ganaste"))
}
ganadores(ejemplo)
