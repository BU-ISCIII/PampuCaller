###############################################################################
# Package: pampuCaller
# Script: misc.R
# 
# Author: smonzon
# date: october 2014
###############################################################################

#' calc_depth function
#'
#' This function calculates depth from variation dataframe. It needs raw read counts in columns named as A,C,G,T for forward reads
#' and a,c,g,t for reverse reads.
#' @param data.frame 
#' @keywords variation
#' @export
#' @examples
#' calc_depth(data.frame)
calc_depth <- function(variation){
	variation$depth <- variation$A + variation$a + variation$C + 
									   + variation$c + variation$T + variation$t + variation$G + variation$g + variation$DEL + variation$INS
	variation
}

#' calc_freq function
#'
#' This function calculates nucleotide frequency. It needs raw read counts in columns named as A,C,G,T for forward reads
#' and a,c,g,t for reverse reads, and a column with position depth called depth.
#' @param data.frame
#' @keywords frequency
#' @export
#' @examples
#' graph_barerror()
calc_freq <- function(variation){
	variation$per_A <- (variation$A + variation$a)/variation$depth
	variation$per_C <- (variation$C + variation$c)/variation$depth
	variation$per_T <- (variation$T + variation$t)/variation$depth
	variation$per_G <- (variation$G + variation$g)/variation$depth
	variation		
}

#' pos_filter function
#'
#' @param data.frame
#' @keywords regions filter
#' @export
#' @examples
#' pos_filter()
pos_filter <- function(variation,start_end,depth){
	
	variation <- variation[variation$POS >= as.numeric(start_end[2]) & variation$POS <= as.numeric(start_end[3]),]	
	variation <- variation[variation$depth >= depth,]
	variation
}