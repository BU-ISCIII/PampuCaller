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
					variation = data.frame(),
					genes = data.frame(),
					polymorphisms = data.frame(),
					samples = data.frame(),
					... ){
						.Object<-callNextMethod(.Object,
						variation = variation ,
						genes = genes,
						polymorphisms = polymorphisms,
						samples = samples,	
						...)				
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

.draw_curve <- function(object,){

}

setMethod("draw_curve",signature="controlSet",.draw_curve)

.mean_sd <- function(object,position,samples,nucleotide){

}

setMethod("mean_sd",signature="controlSet",.mean_sd)