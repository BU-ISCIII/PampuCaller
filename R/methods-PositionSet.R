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
					if(nrow(control) != 0 ){
						control <- calc_depth(control)
						control <- calc_freq(control)
					}

					if(nrow(test) != 0){
						test <- calc_depth(test)
						test <- calc_freq(test)
					}

					#Initialize object
					.Object@control <- control
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
		object <- regions_filter(object,depth=depth)
		control_mean <- mean_sd(object)
		compare_control_set <- new("ControlCompareSet",meanControl=control_mean,test=object@test,samples=object@samples,regions=object@regions,polymorphisms=object@polymorphisms)
	return(compare_control_set)

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
.mean_sd <- function(object,position=NULL,samples=NULL,nucleotide=NULL){
	data <- object@control

	control_mean <- ddply(data,.(POS),function(x){
			ref <- x$REF[1]
			num_samples<-length(x$sample)
			depth<-mean(x$depth)
			mean_per_A<-mean(x$per_A)
			sd_per_A<-sd(x$per_A)
			cv_per_A<- sd_per_A/mean_per_A
			mean_per_C<-mean(x$per_C)
			sd_per_C<-sd(x$per_C)
			cv_per_C<- sd_per_C/mean_per_C
			mean_per_T<-mean(x$per_T)
			sd_per_T<-sd(x$per_T)
			cv_per_T<- sd_per_T/mean_per_T
			mean_per_G<-mean(x$per_G)
			sd_per_G<-sd(x$per_G)
			cv_per_G<- sd_per_G/mean_per_G
			mean_per_a<-mean(x$per_a)
			sd_per_a<-sd(x$per_a)
			cv_per_a<- sd_per_a/mean_per_a
			mean_per_c<-mean(x$per_c)
			sd_per_c<-sd(x$per_c)
			cv_per_c<- sd_per_c/mean_per_c
			mean_per_t<-mean(x$per_t)
			sd_per_t<-sd(x$per_t)
			cv_per_t<- sd_per_t/mean_per_t
			mean_per_g<-mean(x$per_g)
			sd_per_g<-sd(x$per_g)
			cv_per_g<- sd_per_g/mean_per_g
			mean_per_total_A<-mean(x$per_total_A)
			sd_per_total_A<-sd(x$per_total_A)
			cv_per_total_A <- sd_per_total_A/mean_per_total_A
			mean_per_total_C<-mean(x$per_total_C)
			sd_per_total_C<-sd(x$per_total_C)
			cv_per_total_C <- sd_per_total_C/mean_per_total_C
			mean_per_total_T<-mean(x$per_total_T)
			sd_per_total_T<-sd(x$per_total_T)
			cv_per_total_T <- sd_per_total_T/mean_per_total_T
			mean_per_total_G<-mean(x$per_total_G)
			sd_per_total_G<-sd(x$per_total_G)
			cv_per_total_G <- sd_per_total_G/mean_per_total_G
			data.frame(ref=ref,num_samples=num_samples,depth=depth,mean_per_A=mean_per_A,sd_per_A=sd_per_A,cv_per_A=cv_per_A,mean_per_C=mean_per_C,sd_per_C=sd_per_C,cv_per_C=cv_per_C,mean_per_T=mean_per_T,sd_per_T=sd_per_T,cv_per_T=cv_per_T,mean_per_G=mean_per_G,sd_per_G=sd_per_G,cv_per_G=cv_per_G,
				mean_per_a=mean_per_a,sd_per_a=sd_per_a,cv_per_a=cv_per_a,mean_per_c=mean_per_c,sd_per_c=sd_per_c,cv_per_c=cv_per_c,mean_per_t=mean_per_t,sd_per_t=sd_per_t,cv_per_t=cv_per_t,mean_per_g=mean_per_g,sd_per_g=sd_per_g,cv_per_g=cv_per_g,
				mean_per_total_A=mean_per_total_A,sd_per_total_A=sd_per_total_A,cv_per_total_A=cv_per_total_A,mean_per_total_C=mean_per_total_C,sd_per_total_C=sd_per_total_C,cv_per_total_C=cv_per_total_C,mean_per_total_T=mean_per_total_T,sd_per_total_T=sd_per_t,cv_per_total_T=cv_per_total_T,mean_per_total_G=mean_per_total_G,sd_per_total_G=sd_per_total_G,cv_per_total_G=cv_per_total_G)
	})
	return(control_mean)
}
setMethod("mean_sd",signature="PositionSet",.mean_sd)

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