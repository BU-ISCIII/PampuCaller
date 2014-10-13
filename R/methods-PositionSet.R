###############################################################################
# Package: pampuCaller
# Script: methods-PositionSet.R
# 
# Author: smonzon
###############################################################################

######
## Initialize 
#####

setMethod("initialize","PositionSet",
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
					sample <- unique(variation$sample)
					.Object@samples <- as.data.frame(sample)
					.Object@variation <- variation
					return(.Object)
			}
)	

#######
## Class Validation
#########
## Implement real validation
setValidity("PositionSet",function(object){
		TRUE
		}
)		

#######
## Class Methods
#######
setMethod("show","PositionSet",
		function(object){
			size<-dim(object)
			cat(class(object), "instance with:\n\t",size[1],"variant positions\n\t",size[2],"genes\n\t",size[3],"polymorphisms\n\t",size[4],"samples\n")
		}
)

setMethod("dim","PositionSet",
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

#' Graph_barerror function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' graph_barerror()
graph_barerror <- function(object,depth,samples=object@samples$sample,zoom,...){
	var <- get_variation(object)
	var <- var[var$depth > depth & var$sample %in% samples,]
	var_ggplot <- melt(var,id.vars=c('sample','POS'), measure.vars=c('per_A','per_G','per_T','per_C'))

	g <-ggplot(var_ggplot, aes(POS,value, fill=variable)) +
 		geom_bar(stat="identity",position="dodge") + 
 		ylim(zoom) +
 		theme_bw() + 
 		theme(axis.text.x = element_text(size = 5.5,angle=75, vjust=0.5), strip.text.x = element_text(size=6.5)) + 
 		scale_fill_manual(values=cbPalette) + 
 		labs(title="Amplicon Position Error", x="POS",y="Alt Freq",legend="Nucleotide percentage") +
 		facet_wrap(~sample)
	g
} 
# 
# scale_y_continuous(trans=log2_trans()) +	
setMethod("graph_barerror",signature="PositionSet",graph_barerror)

#' mean_sd function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' mean_sd()
mean_sd <- function(object,position,samples,nucleotide){

}

setMethod("mean_sd",signature="PositionSet",mean_sd)

#' get_variation function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' get_variation()
get_variation <- function(object){
	variation <- object@variation
	variation
}
setMethod("get_variation",signature="PositionSet",get_variation)

#' get_sample function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' get_sample()
get_sample <- function(object){
	samples <- object@samples
	samples
}
setMethod("get_sample",signature="PositionSet",get_sample)