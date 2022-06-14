library(BoolNet)
#generar las redes
sink("activador_con_glia")
cat("targets, factors\n")
cat("Pre, Pre \n")
cat("Pos, Pre \n")
cat("Glia, Pre \n")
sink()

sink("inhibidor_con_glia")
cat("targets, factors\n")
cat("Pre, Pre \n")
cat("Pos, ! Pre \n")
cat("Glia, Pre \n")
sink()

sink("modelo1_3_celulas")
cat("targets, factors\n")
cat("Cel1, Cel1 \n")
cat("Cel2, Cel1 \n")
cat("Cel3, Cel2 \n")
sink()

sink("modelo2_3_celulas")
cat("targets, factors\n")
cat("Cel1, Cel1 \n")
cat("Cel2, ! Cel1 \n")
cat("Cel3, Cel2 \n")
sink()

sink("modelo3_3_celulas")
cat("targets, factors\n")
cat("Cel1, Cel1 \n")
cat("Cel2, Cel1 \n")
cat("Cel3, ! Cel2 \n")
sink()

sink("modelo4_3_celulas")
cat("targets, factors\n")
cat("Cel1, Cel1 \n")
cat("Cel2, ! Cel1 \n")
cat("Cel3, ! Cel2 \n")
sink()

#cargar las redes en objetos
activador.glia<-loadNetwork("activador_con_glia")
inhibidor.glia<-loadNetwork("inhibidor_con_glia")
modelo_1<-loadNetwork("modelo1_3_celulas")
modelo_2<-loadNetwork("modelo2_3_celulas")
modelo_3<-loadNetwork("modelo3_3_celulas")
modelo_4<-loadNetwork("modelo4_3_celulas")

#obtener datos de atractores y posibles estados de las redes
getAttractors(activador.glia)
plotAttractors(getAttractors(activador.glia))
plotStateGraph(getAttractors(activador.glia),drawLabels = T)

getAttractors(inhibidor.glia)
plotAttractors(getAttractors(inhibidor.glia))
plotStateGraph(getAttractors(inhibidor.glia),drawLabels = T)

getAttractors(modelo_1)
plotAttractors(getAttractors(modelo_1))
plotStateGraph(getAttractors(modelo_1),drawLabels = T)

getAttractors(modelo_2)
plotAttractors(getAttractors(modelo_2))
plotStateGraph(getAttractors(modelo_2),drawLabels = T)

getAttractors(modelo_3)
plotAttractors(getAttractors(modelo_3))
plotStateGraph(getAttractors(modelo_3),drawLabels = T)

getAttractors(modelo_4)
plotAttractors(getAttractors(modelo_4))
plotStateGraph(getAttractors(modelo_4),drawLabels = T)



