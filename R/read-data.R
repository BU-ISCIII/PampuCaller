###############################################################################
# Package: pampuCaller
# Script: read-data.R
# 
# Author: smonzon
# date: October 2014
###############################################################################

#' read_targets() function
#'
#' This function reads targets file
#' @param targets file
#' @keywords targets
#' @export
#' @examples
#' read_targets(file="targets.txt")
read_targets <- function(file="targets.txt"){
	t <- read.table(file,sep="\t", header=T)
	####
	## TO DO: Format control
	###
	return(t)
}

#' read_regions() function
#'
#' This function reads regions file
#' @param regions file
#' @keywords regions
#' @export
#' @examples
#' read_targets(file="regions.bed")
read_regions <- function(file="regions.bed"){
	r <- read.table(file,sep="\t")
	colnames(r) <- c("chr","start","end","name","strand")
	####
	## TO DO: Format control
	###
	return(r)
}

#' create_position_set() function
#'
#' This function uses targets, regions and polymorphism info and creates a PositionSet instance
#' @param regions bed file
#' @keywords PositionSet instance
#' @export
#' @examples
#' create_position_set(targets,regions,polymorphisms)
create_position_set <- function(targets=data.frame(),regions=data.frame(),polymorphisms=data.frame()){

	control <- apply(as.matrix(targets[targets$ControlTest == "CONTROL",]),1,read_data)
	test <- apply(as.matrix(targets[targets$ControlTest == "TEST",]),1,read_data)
	control <- do.call(rbind,control)
	test <- do.call(rbind,test)

	position_set <- new("PositionSet",control=control,test=test,samples=targets,regions=regions,polymorphisms=polymorphisms)
	return(position_set)
}

read_data <- function(targ){
	f <- read.table(targ["filename"],header=T,sep="\t")
	f$sample <- targ["Sample"]
	return(f)
}