###############################################################################
# Package: pampuCaller
# Script: misc.R
# 
# Author: smonzon
# date: october 2014
###############################################################################

## Graphs palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

brewer_qualitative <- c("#0000ff","#ff0000","#483215","#008900","#7244c4","#e65a11","#000000","#e6e528","#ff00ee","#6e0000","#00c7dd","#d4b455","#8f008d","#736b00","#7d8cbf")


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
									   + variation$c + variation$T + variation$t + variation$G + variation$g
	
	## Calc all nt totals
	variation$total_A <- (variation$A + variation$a)
	variation$total_C <- (variation$C + variation$c)
	variation$total_T <- (variation$T + variation$t)
	variation$total_G <- (variation$G + variation$g)

	## TODO Not counting indels forward and reverse separated.
	variation$forward_reads <- variation$A + variation$C + variation$T + variation$G 
	variation$reverse_reads <- variation$a + variation$c + variation$t + variation$g
	variation$per_forward_reads <- variation$forward_reads/variation$depth
	variation$per_reverse_reads <- variation$reverse_reads/variation$depth
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
	# Forward freq per nt
	variation$per_A <- (variation$A)/variation$depth
	variation$per_C <- (variation$C)/variation$depth
	variation$per_T <- (variation$T)/variation$depth
	variation$per_G <- (variation$G)/variation$depth
	# Reverse freq per nt
	variation$per_a <- (variation$a)/variation$depth
	variation$per_c <- (variation$c)/variation$depth
	variation$per_t <- (variation$t)/variation$depth
	variation$per_g <- (variation$g)/variation$depth

	# Total freq per nt
	variation$per_total_A <- (variation$A + variation$a)/variation$depth
	variation$per_total_C <- (variation$C + variation$c)/variation$depth
	variation$per_total_T <- (variation$T + variation$t)/variation$depth
	variation$per_total_G <- (variation$G + variation$g)/variation$depth
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
	variation <- variation[variation$depth >= as.numeric(depth),]
	variation
}