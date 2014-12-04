cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

brewer_qualitative <- c("#0000ff","#ff0000","#483215","#008900","#7244c4","#e65a11","#000000","#e6e528","#ff00ee","#6e0000","#00c7dd","#d4b455","#8f008d","#736b00","#7d8cbf")

## Test una muestra con dos réplicas

rep1 <- read.table("data/IMR90_EX14_REP1-align.bam-count_variants.vcf",header=T)
rep2 <- read.table("data/IMR90_EX14_REP2-align.bam-count_variants.vcf",header=T)
rep1$sample <- "IMR90_EX14_REP1"                                                
rep2$sample <- "IMR90_EX14_REP2"                                               
variants <- rbind(rep1,rep2)
variants_o <- new("PositionSet",variation=variants)

graph_barerror(variants_o,10000)

#####################################################################

## Test directorio con muchas muestras.
create_dataframe <- function(filename){
	names <- strsplit(filename,"/")
	reads_file_name <- names[[1]][length(names[[1]])]
	sample <- strsplit(reads_file_name,"-")[[1]][1]

	var <- read.table(filename,sep="\t",header=T)
	var$sample <- sample
	return(var)
}

input_dir <- "/home/smonzon/Documentos/proyectos/20140304_RBMOSAIC/ANALYSIS/20140327_RBMOSAIC02_DEEPSEQ/count_variants"
result_dir <- "/home/smonzon/Documentos/proyectos/20140304_RBMOSAIC/ANALYSIS/20140327_RBMOSAIC02_DEEPSEQ/graphics_R"
pattern <- "count_variants.vcf"

reads_file <- list.files(input_dir,recursive=TRUE,pattern=pattern,full.names=TRUE)
resul <- lapply(reads_file,create_dataframe)
resul_df <- do.call(rbind,resul)
resul_df <- as.data.frame(resul_df)

variants_o <- new("controlSet",variation=resul_df)

samples <- get_sample(variants_o)

recorrer_muestras(samples,"IMR90",result_dir)


recorrer_muestras <- function(samples,pattern_control,dir){
	control <- samples[grep("IMR90.*REP1",samples$sample,perl=T),]
	control_m <- do.call(rbind,strsplit(as.character(control),"_"))

	apply(control_m,1,draw_graph,samples=samples,pattern=pattern_control,dir)
}

draw_graph <- function(exon,samples,pattern_control,dir){
		print(paste("Graph",exon[2]))
		name_c <- samples[grep(paste(pattern_control,".*",exon[2],sep=""),samples$sample,perl=T),]
		name <- samples[grep(paste("[^",pattern_control,"]_",exon[2],sep=""),samples$sample,perl=T),]

		dir <- paste(dir,"/",exon[2],sep="")

		dir.create(dir)

		pdf(paste(dir,"/",exon[2],"_zoom_0.015",sep=""),width=40)
		print(graph_barerror(variants_o,10000,name_c,zoom=c(0,0.015)))
	    print(graph_barerror(variants_o,10000,name,zoom=c(0,0.015)))
		dev.off()

		pdf(paste(dir,"/",exon[2],"_zoom_0.008",sep=""),width=40)
		print(graph_barerror(variants_o,10000,name_c,zoom=c(0,0.008)))
	    print(graph_barerror(variants_o,10000,name,zoom=c(0,0.008)))
		dev.off()

		pdf(paste(dir,"/",exon[2],"_zoom_0.5",sep=""),width=40)
		print(graph_barerror(variants_o,10000,name_c,zoom=c(0,0.5)))
	    print(graph_barerror(variants_o,10000,name,zoom=c(0,0.5)))
		dev.off()

		pdf(paste(dir,"/",exon[2],"zoom_0.25",sep=""),width=40)
		print(graph_barerror(variants_o,10000,name_c,zoom=c(0,0.25)))
	    print(graph_barerror(variants_o,10000,name,zoom=c(0,0.25)))
		dev.off()
	}


## Test gráficas para 4 controles en 454 que tienen los 27 exones de RB.

input_dir <- "/home/smonzon/Documentos/proyectos/20140304_RBMOSAIC/ANALYSIS/20140327_RBMOSAIC03_DEEPSEQ/count_variants"
result_dir <- "/home/smonzon/Documentos/proyectos/20140304_RBMOSAIC/ANALYSIS/20140327_RBMOSAIC03_DEEPSEQ/graphics_R"
pattern <- "count_variants.vcf"

reads_file <- list.files(input_dir,recursive=TRUE,pattern=pattern,full.names=TRUE)
resul <- lapply(reads_file,create_dataframe)
resul_df <- do.call(rbind,resul)
resul_df <- as.data.frame(resul_df)

variants_o <- new("controlSet",variation=resul_df)

# Graphs
pdf("foursamplescontrol_zoom_0.015",width=60,height=30)
print(graph_barerror(variants_o,1000,zoom=c(0,0.015)))
dev.off()

pdf("foursamplescontrol_zoom_0.008",width=60,height=30)
print(graph_barerror(variants_o,1000,zoom=c(0,0.008)))
dev.off()

pdf("foursamplescontrol_zoom_0.5",width=60,height=30)
print(graph_barerror(variants_o,1000,zoom=c(0,0.5)))
dev.off()

pdf("foursamplescontrol_zoom_0.25",width=60,height=30)
print(graph_barerror(variants_o,1000,zoom=c(0,0.25)))
dev.off()


################################################
### PRUEBA Octubre 2014
############ ####################################

library(pampuCaller)
targets <- read_targets()
regions <- read_regions("amplicones.bed")
position_set <- create_position_set(targets,regions) 
compare_control_set <- calling_prep(position_set,depth=1000)


## TODO
recorrer_amplicones <- function(reg,dir){
	control <- samples[grep("IMR90.*REP1",samples$sample,perl=T),]
	control_m <- do.call(rbind,strsplit(as.character(control),"_"))

	apply(control_m,1,draw_graph,samples=samples,pattern=pattern_control,dir)
}

draw_graph <- function(exon,samples,pattern_control,dir){
		print(paste("Graph",exon[2]))
		name_c <- samples[grep(paste(pattern_control,".*",exon[2],sep=""),samples$sample,perl=T),]
		name <- samples[grep(paste("[^",pattern_control,"]_",exon[2],sep=""),samples$sample,perl=T),]

		dir <- paste(dir,"/",exon[2],sep="")

		dir.create(dir)

		pdf(paste(dir,"/",exon[2],"_zoom_0.015",sep=""),width=40)
		print(graph_barerror(variants_o,10000,name_c,zoom=c(0,0.015)))
	    print(graph_barerror(variants_o,10000,name,zoom=c(0,0.015)))
		dev.off()

		pdf(paste(dir,"/",exon[2],"_zoom_0.008",sep=""),width=40)
		print(graph_barerror(variants_o,10000,name_c,zoom=c(0,0.008)))
	    print(graph_barerror(variants_o,10000,name,zoom=c(0,0.008)))
		dev.off()

		pdf(paste(dir,"/",exon[2],"_zoom_0.5",sep=""),width=40)
		print(graph_barerror(variants_o,10000,name_c,zoom=c(0,0.5)))
	    print(graph_barerror(variants_o,10000,name,zoom=c(0,0.5)))
		dev.off()

		pdf(paste(dir,"/",exon[2],"zoom_0.25",sep=""),width=40)
		print(graph_barerror(variants_o,10000,name_c,zoom=c(0,0.25)))
	    print(graph_barerror(variants_o,10000,name,zoom=c(0,0.25)))
		dev.off()
	}



## Test
position_set_filt <- regions_filter(position_set,depth=1000)
pfc <- get_control(position_set_filt)
a <- mean_sd(position_set_filt)

start_end <- get_regions(compare_control_set)[12,]
depth <- 0
zoom <- c(0,0.008)
samples <- compare_control_set@samples$Sample
var <- get_test(compare_control_set)
var_c <- get_control(compare_control_set)

l <- graph_barerror(compare_control_set,start_end=start_end,zoom=c(0,0.008))
