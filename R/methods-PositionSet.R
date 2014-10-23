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

					#Calc precomputations
					control <- calc_depth(control)
					control <- calc_freq(control)

					test <- calc_depth(test)
					test <- calc_freq(test)

					#Initialize object
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
## TO DO: Implement real validation
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

#' calling_prep function
#' This function calculates 
#' 
#' @param object
#' @param depth
#' @param calling c(control,filter)
#' @keywords region filter
#' @export
#' @examples
#' region_filter(object)
.calling_prep <- function(object,depth=20,calling="control"){
	if(calling == "control"){
		object <- regions_filter(object)
		object <- mean_sd(object)

	}else if(calling == "filter"){

	}
}
setMethod("calling_prep",signature="PositionSet",.calling_prep)


#' regions_filter function
#' This function filters positions in control and test slot from positionSet per regions in regions slot
#' 
#' @param object
#' @keywords region filter
#' @export
#' @examples
#' region_filter(object)
.regions_filter <- function(object,depth=20){
	control_fil <- apply(as.matrix(object@regions),1,pos_filter,variation=object@control,depth=depth)
	test_fil <- apply(as.matrix(object@regions),1,pos_filter,variation=object@test,depth=depth)
	control_fil <- do.call(rbind,control_fil)
	test_fil <- do.call(rbind,test_fil)
	object@control <- control_fil
	object@test <- test_fil
	return(object)
}
setMethod("regions_filter",signature="PositionSet",.regions_filter)


#' mean_sd function
#' This function creates a ControlCompareSet when control slot has info.	
#' 
#' @param object
#' @param samples
#' @keywords mean
#' @export
#' @examples
#' mean_sd(object)
.mean_sd <- function(object,position,samples,nucleotide){

}
setMethod("mean_sd",signature="PositionSet",.mean_sd)


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
.get_control <- function(object){
	v <- object@control
	v
}
setMethod("get_control",signature="PositionSet",.get_control)

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
setMethod("get_test",signature="PositionSet",.get_test)

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
setMethod("get_samples",signature="PositionSet",.get_samples)

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
setMethod("get_regions",signature="PositionSet",.get_regions)