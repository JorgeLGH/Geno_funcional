#----Ejercicio 1----
tab_1 <- read.table("expr_normalizada (1).txt", header = T, as.is = T)#cargar la tabla

boxplot(tab_1)#boxplot y se ve normalizado

library(Biobase)
library(limma)
pca <- prcomp(tab_1)#generar el pca
biplot(pca) #plotear el pca, tarda pero se hace

expe <- tab_1[which(colnames(tab_1)==colnames(tab_1)[c(1:3)])]
control<- tab_1[which(colnames(tab_1)==colnames(tab_1)[c(4:6)])]

types <- factor(c("EX","EX","EX","WT","WT","WT"))
types
design <- model.matrix(~ 0+types)
colnames(design) <- levels(types)
design
contMatrix <- makeContrasts(EX-WT, levels=design)
contMatrix
fit  <- lmFit(tab_1,design)
fit2 <- contrasts.fit(fit,contMatrix)
fit2 <- eBayes(fit2)
deTable <- topTable(fit2, number=nrow(fit2), lfc=log2(1.5), p.value=0.05)#aquí ya estoy poniendo mis criterios para decir si están diferencialmente expresados
#o no, basado en significancia estadística y cambio en magnitud de expresión (log-fold change), se puede hacer con todos si le quito los parámetros que delimitan
#cuáles datos entran y cuáles no
sum(deTable$logFC>0)#los de mayor expresión
sum(deTable$logFC<0)#los de menor expresión


fullTable <- topTable(fit2, number=nrow(fit2), sort.by="logFC")
plot(x=fullTable$logFC, y=fullTable$B)#primera vez que me sale esta cosa, creo, pero no estoy así como que muy seguro porque no me hace sentido 
#el -5 en el eje Y pero pues es lo que hay
#claramente hay una cantidad enorme de genes, de los cuales parece ser que hay relativameente un número bajo de ellos que se encuentra diferencialmente
#expresados, lo cual tiene sentido. También tiene sentido que en su mayoría no estén diferencialmente expresados, aunque no conozco el experimento para
#poder decirlo con seguridad

#----Ejercicio 2----
#Breve introducción a IRanges
library(Biobase)
library(IRanges)
x <- IRanges(start=c(11,35,40), end=c(20,50,63))#crar los rangos, inicio y fin

#datos generales que pueden ser usados
start(x) # Todos los inicios
end(x)   # Todos los finales
width(x) # Ancho de cada rango
range(x) # Rango total de todos

coverage(x) # La suma de la cobertura en cada posición
reduce(x) # une los rangos encimados, no aparecen los que son "repetidos"

exons <- reduce(x)
reads <- IRanges(start=c(1,21,30,50,80), width=20)
reads
countOverlaps(exons, reads)#encontrar donde se sobrelapan las lecturas de mis reads y los exones

#Obteniendo la anotación de un genoma
load("human.Rdata")#descargué este para más rápido, datos de anotación de humano
human

#info básica que podemos sacar, todas dicen lo que el nombre indica
seqnames(human)
ranges(human)
strand(human)
mcols(human)

table(mcols(human)$gene_biotype)#tabla resumiendo información
mcols(human) <- mcols(human)[,c("source","gene_id","gene_name","gene_biotype")]
mcols(human)#más fácil manejo, dataframe con los datos organizados

#1.-¿Cómo le harían para quedarse exclusivamente con las anotaciones de "miRNA"?
(human)[which(human$source=="miRNA")]
#2.-¿y solamente aquellas anotaciones de la cadena "-"?
human[which(strand(human)=="-")]

#Anotación de secuencias mapeadas
library(Rsamtools)

what <- c("rname", "strand", "pos", "qwidth")
param <- ScanBamParam(what=what)
bam <- scanBam("human_mapped.bam", param=param)#con los parámetros de arriba, solo extraemos lo que queremos
class(bam)
lapply(bam, names)

mapGR <- GRanges(
  seqnames = bam[[1]]$rname,
  ranges   = IRanges(start=bam[[1]]$pos, width=bam[[1]]$qwidth),
  strand   = bam[[1]]$strand
)
mapGR #con estoy ya tenemos los GRanges de los datos que teníamos
mcols(human)$counts = countOverlaps(human, mapGR) #es una columna que usa el countoverlaps, aunque los rangos podrían estar
#contenidos unos dentro de otros, lo que elimina info
mcols(human)
typeCounts <- aggregate(mcols(human)$counts, by=list("biotype"=mcols(human)$gene_biotype), sum)
typeCounts #el conteo de los objetos dada la categoría que querramos
geneCounts <- aggregate(mcols(human)$counts, by=list("id"=mcols(human)$gene_name), sum)
head(geneCounts) #lo mismo que el de arriba pero con genes individuales

minCount <- 40000
typeCountsHigh <- typeCounts[typeCounts$x > minCount,]
typeCountsHigh <- typeCountsHigh[order(typeCountsHigh$x),]
typeCountsHigh <- rbind(data.frame("biotype"="other",
                                  "x"=sum(typeCounts$x[typeCounts$x <= minCount])),
                       typeCountsHigh)

pie(typeCountsHigh$x, labels=typeCountsHigh$biotype, col=rev(rainbow(nrow(typeCountsHigh))),
    main="Number of aligned reads per biotype")#visualización de los distintos tipos de "biotype" de las secuencias analizadas 

#----Ejercicio 3----
library(sleuth)
tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}
t2g <- tx2gene()

base_dir<-"Archivo/"

samples<- paste0("sample", c(2:12))
muestra<- paste0("sample", c(1,1,1,1,1,1,1,1,1,2,2))#creo que puedo saber cuáles son de qué experimento repitiendo el análisis
#aunque sean diferentes muestras, los resultados deben de señalar un parecido alto en el análisis si corresponden al mismo grupo, CREO; si no es así
#pues ya no tengo tiempo para probar otra cosa ggg (es que tarda mucho). Mi teoría es que la parte interactiva me dará las muestras gráficas de la similitud
#entre las muestras
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, muestras = muestra, stringsAsFactors=FALSE)
so <- sleuth_prep(s2c, ~sample, target_mapping = t2g, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~muestras, 'full')
so <- sleuth_wt(so, which_beta="muestrassample1")                   
sleuth_live(so)
#si mi expermiento sirvió y lo interpreté de manera correcta, los que son de un mismo grupo, son del 1:6
#creo que iba por buen camino, pero se me agotó el tiempo :(

samples<- paste0("sample", c(1:4))
muestra<- paste0("sample", c(8,9,10,10))
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, muestras = muestra, stringsAsFactors=FALSE)
so <- sleuth_prep(s2c, ~sample, target_mapping = t2g, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~muestras, 'full')
so <- sleuth_wt(so, which_beta="muestrassample8")                   
sleuth_live(so)
#no pues creo que está mal hecho jajaj
setwd('C:/Users/fotgo/OneDrive/Documentos/R/Genomica_funcional/')
resultados<-read.table("test_table.csv",sep=",",
                       header=TRUE)
significativos<-which(resultados$qval<0.1)
significativos<-resultados[significativos,]
upregulated<-which(significativos$b>0)
upregulated<-significativos[upregulated,]
downregulated<-which(significativos$b<0)
downregulated<-significativos[downregulated,]
write.table(upregulated,file="C:/Users/fotgo/OneDrive/Documentos/R/Genomica_funcional/Upregulated.txt",sep="\t")
write.table(downregulated,file="C:/Users/fotgo/OneDrive/Documentos/R/Genomica_funcional/Downregulated.txt",sep="\t")
#se supone que sale así con GO en panther, pero la verdad ya no estoy seguro de nada de lo que hice
#1	cellular process (GO:0009987)	51	100.0%	42.9%
#2	metabolic process (GO:0008152)	34	66.7%	28.6%
#3	biological regulation (GO:0065007)	17	33.3%	14.3%
#4	localization (GO:0051179)	17	33.3%	14.3%
#----Ejercicio 4----
library(BoolNet)
sink("red_ex")
cat("targets, factors\n")
cat("B, ! A\n")
cat("A, ! C\n")
cat("C, ! B\n")
sink()#crear la red
red_ex <- loadNetwork("red_ex")
getAttractors(red_ex) #hay 2 atractores
#el estado más probable depende del estado inicial, pero hay más posibilidades que se esté oscilando en los estados
#asociados al atractor cíclico debido a que tienen una cuenca de atractor mayor
plotAttractors(getAttractors(red_ex))#representación de atractores
plotStateGraph(getAttractors(red_ex))#representación de todos los atractores
