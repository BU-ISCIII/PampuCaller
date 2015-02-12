###############################################################################
# Package: pampuCaller
# Script: methods-VariationSet.R
# 
# Author: smonzon
# date: october 2014
###############################################################################

######
## Initialize 
#####
setMethod("initialize","VariationSet",
			function(.Object,
					variants = data.frame(),
					annotation = data.frame(),
					... ){		

					#Initialize object
					.Object@variants <- variants
					.Object@annotation <- annotation
					return(.Object)
			}
)	

#######
## Class Validation
#########
## TO DO: Implement real validation
setValidity("VariationSet",function(object){
		TRUE
		}
)		

#######
## Class Methods
#######
setMethod("show","VariationSet",
		function(object){
			size<-dim(object)
			cat(class(object), "instance with:\n\t",size[1],"Variant positions\n\t")
		}
)

setMethod("dim","VariantSet",
		function(x){
			size_variants <- nrow(x@variants)
		}
)

############
### Accession methods
############


#' get_variants function
#'
#' Accesion function to control dataframe
#' @param object
#' @keywords control variation
#' @export
#' @examples
#' get_test(object)
.get_variants <- function(object){
	v <- object@variants
	v
}
setMethod("get_variants",signature="VariantSet",.get_variants)