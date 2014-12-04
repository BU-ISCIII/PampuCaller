###############################################################################
# Package: pampuCaller
# Script: methods-ControlCompareSet.R
# 
# Author: smonzon
# date: october 2014
###############################################################################

######
## Initialize 
#####

setMethod("initialize","ControlCompareSet",
			function(.Object,
					meanControl = data.frame(),
					test = data.frame(),
					regions = data.frame(),
					polymorphisms = data.frame(),
					samples = data.frame(),
					... ){		

					#Initialize object
					.Object@meanControl <- meanControl
					.Object@test <- test
					.Object@regions <- regions
					.Object@polymorphisms <- polymorphisms
					.Object@samples <- samples
					return(.Object)
			}
)	

#######
## Class Validation
#########
## TO DO: Implement real validation
setValidity("ControlCompareSet",function(object){
		TRUE
		}
)		

#######
## Class Methods
#######
setMethod("show","ControlCompareSet",
		function(object){
			size<-dim(object)
			cat(class(object), "instance with:\n\t",size[1],"Mean control positions\n\t",size[2],"test positions\n\t",size[3],"regions\n\t",size[4],"polymorphisms\n\t",size[5],"samples\n")
		}
)

setMethod("dim","ControlCompareSet",
		function(x){
			size_control <- nrow(x@meanControl)
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
#' @keywords barerror
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


############
### Accession methods
############


#' get_control function
#'
#' Accesion function to control dataframe
#' @param object
#' @keywords control variation
#' @export
#' @examples
#' get_test(object)
.get_mean_control <- function(object){
	v <- object@meanControl
	v
}
setMethod("get_control",signature="ControlCompareSet",.get_mean_control)

#' get_test function
#'
#' Accesion function to test dataframe
#' @param object
#' @keywords test variation
#' @export
#' @examples
#' get_variation(object)
.get_test <- function(object){
	v <- object@test
	v
}
setMethod("get_test",signature="ControlCompareSet",.get_test)

#' get_samples function
#'
#' Accession function for samples
#' @keywords sample
#' @export
#' @examples
#' get_samples(object)
.get_samples <- function(object){
	s <- object@samples
	s
}
setMethod("get_samples",signature="ControlCompareSet",.get_samples)

#' get_regions function
#'
#' Accession function for regions
#' @keywords regions
#' @export
#' @examples
#' get_regions(object)
.get_regions <- function(object){
	s <- object@regions
	s
}
setMethod("get_regions",signature="ControlCompareSet",.get_regions)