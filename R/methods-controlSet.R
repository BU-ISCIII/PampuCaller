###############################################################################
# Package: pampuCaller
# Script: methods-controlSet.R
# 
# Author: smonzon
###############################################################################

######
## Initialize 
#####

setMethod("initialize","controlSet",
			function(.Object,
					variation = .Object@variation,
					genes = data.frame(),
					polymorphisms = data.frame(),
					samples = data.frame(),
					... ){		

					variation$depth <- variation$A + variation$a + variation$C + 
									   + variation$c + variation$T + variation$t + variation$G + variation$g + variation$DEL + variation$INS
					variation$per_A <- (variation$A + variation$a)/variation$depth
					variation$per_C <- (variation$C + variation$c)/variation$depth
					variation$per_T <- (variation$T + variation$t)/variation$depth
					variation$per_G <- (variation$G + variation$g)/variation$depth
					.Object@variation <- variation
					return(.Object)
			}
)	

#######
## Class Validation
#########
## Implement real validation
setValidity("controlSet",function(object){
		TRUE
		}
)		

#######
## Class Methods
#######
setMethod("show","controlSet",
		function(object){
			size<-dim(object)
			cat(class(object), "instance with:\n\t",size[1],"variant positions\n\t",size[2],"genes\n\t",size[3],"polymorphisms\n\t",size[4],"samples\n")
		}
)

setMethod("dim","controlSet",
		function(x){
			size_variants <- nrow(x@variation)
			size_genes <- nrow(x@genes)
			size_po <- nrow(x@polymorphisms)
			size_samples <- nrow(x@samples)
			c(size_variants,size_genes,size_po,size_samples)
		}
)

##################
## Features
##################

.graph_barerror <- function(object,depth,...){
	var <- get_variation(object)
	var <- var[var$depth > depth,]
	var_ggplot <- melt(var,id.vars=c('sample','POS'), measure.vars=c('per_A','per_G','per_T','per_C'))

	g <-ggplot(var_ggplot, aes(POS,value, fill=variable),scale) +
 		geom_bar(stat="identity",position="dodge") + 		
 		ylim(0,0.0008) +
 		theme(axis.text.x = element_text(size = 5.5,angle=75, vjust=0.5), strip.text.x = element_text(size=6.5)) + 
 		scale_fill_manual(values=cbPalette) + 
 		theme_bw() + 
 		labs(title="Amplicon Position Error", x="POS") +
 		facet_wrap(~sample)
	g
}
#scale_y_continuous(trans=log_trans(),breaks=c(0.01,0.1,0.5,0.9,1.0)) +
setMethod("graph_barerror",signature="controlSet",.graph_barerror)

.mean_sd <- function(object,position,samples,nucleotide){

}

setMethod("mean_sd",signature="controlSet",.mean_sd)

.get_variation <- function(object){
	variation <- object@variation
	variation
}
setMethod("get_variation",signature="controlSet",.get_variation)
