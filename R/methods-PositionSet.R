###############################################################################
# Package: pampuCaller
# Script: methods-PositionSet.R
# 
# Author: smonzon
# date: september 2014
###############################################################################

######
## Initialize 
#####

setMethod("initialize","PositionSet",
			function(.Object,
					control = data.frame(),
					test = data.frame(),
					regions = data.frame(),
					polymorphisms = data.frame(),
					samples = data.frame(),
					... ){		

					control <- calc_depth(control)
					control <- calc_freq(control)

					test <- calc_depth(test)
					test <- calc_freq(test)

					.Object@control <- control
					.Object@test <- test
					.Object@regions <- regions
					.Object@polymorphisms <- polymorphisms
					.Object@samples <- samples
					return(.Object)
			}
)	

calc_depth <- function(variation){
	variation$depth <- variation$A + variation$a + variation$C + 
									   + variation$c + variation$T + variation$t + variation$G + variation$g + variation$DEL + variation$INS
	variation
}

calc_freq <- function(variation){
	variation$per_A <- (variation$A + variation$a)/variation$depth
	variation$per_C <- (variation$C + variation$c)/variation$depth
	variation$per_T <- (variation$T + variation$t)/variation$depth
	variation$per_G <- (variation$G + variation$g)/variation$depth
	variation		
}

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
			cat(class(object), "instance with:\n\t",size[1],"control positions\n\t",size[2],"test positions\n\t",size[3],"regions\n\t",size[4],"polymorphisms\n\t",size[5],"samples\n")
		}
)

setMethod("dim","PositionSet",
		function(x){
			size_control <- nrow(x@control)
			size_test <- nrow(x@test)
			size_regions <- nrow(x@regions)
			size_po <- nrow(x@polymorphisms)
			size_samples <- nrow(x@samples)
			c(size_control,size_test,size_regions,size_po,size_samples)
		}
)

##################
## Features
##################

#' Graph_barerror function
#'
#' This function allows you to draw a a bar graph showing frequency error per position.
#' @param PositionSet Object.
#' @param depth filter
#' @param samples to print
#' @param frequency zoom
#' @keywords 
#' @export
#' @examples
#' graph_barerror()
.graph_barerror <- function(object,depth,samples=object@samples$sample,zoom,...){
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
setMethod("graph_barerror",signature="PositionSet",.graph_barerror)

#' mean_sd function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' mean_sd()
.mean_sd <- function(object,position,samples,nucleotide){

}

setMethod("mean_sd",signature="PositionSet",.mean_sd)

#' get_variation function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' get_variation()
.get_variation <- function(object){
	v <- object@variation
	v
}
setMethod("get_variation",signature="PositionSet",.get_variation)

#' get_sample function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' get_sample()
.get_samples <- function(object){
	s <- object@samples
	s
}
setMethod("get_samples",signature="PositionSet",.get_samples)