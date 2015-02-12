###############################################################################
# Package: pampuCaller
# Script: pampu_caller.R
# 
# Author: smonzon
# date: January 2015
###############################################################################

#' pampu caller function
#'
#' @param 
#' @keywords 
#' @export
#' @examples
#'

pampu_caller_test_1 <- function(x,control,min_freq,times){
	variants <- data.frame()

	seek_control_pos <- control[control["POS"] == as.numeric(x["POS"]),]

	if(is.na(seek_control_pos["mean_per_total_A"]) || is.na(seek_control_pos["mean_per_total_C"]) || is.na(seek_control_pos["mean_per_total_G"]) || is.na(seek_control_pos["mean_per_total_T"])){ 
		print(paste("POS:",x["POS"],"not found in control."))
		return(variants)
	}

	## Test each nt independently 
	if(is.zero(as.numeric(seek_control_pos["mean_per_total_A"]))){
		if(x["per_total_A"] > times/x["depth"]){
			x["ALT"] <- "A"
			variants <- rbind(variants,x)
		}
	}else{
		if(x["per_total_A"] > times*as.numeric(seek_control_pos["mean_per_total_A"])){
			x["ALT"] <- "A"
			variants <- rbind(variants,x)
		}
	}
	if(is.zero(as.numeric(seek_control_pos["mean_per_total_C"]))){
		if(x["per_total_C"] > times/x["depth"]){
			x["ALT"] <- "C"
			variants <- rbind(variants,x)
		}
	}else{
		if(x["per_total_C"] > times*as.numeric(seek_control_pos["mean_per_total_C"])){
			x["ALT"] <- "C"
			variants <- rbind(variants,x)
		}
	}
	if(is.zero(as.numeric(seek_control_pos["mean_per_total_T"]))){
		if(x["per_total_T"] > times/x["depth"]){
			x["ALT"] <- "T"
			variants <- rbind(variants,x)
		}
	}else{
		if(x["per_total_T"] > times*as.numeric(seek_control_pos["mean_per_total_T"])){
			x["ALT"] <- "T"
			variants <- rbind(variants,x)
		}
	}
	if(is.zero(as.numeric(seek_control_pos["mean_per_total_G"]))){
		if(x["per_total_G"] > times/x["depth"]){
			x["ALT"] <- "G"
			variants <- rbind(variants,x)
		}
	}else{
		if(x["per_total_G"] > times*as.numeric(seek_control_pos["mean_per_total_G"])){
			x["ALT"] <- "G"
			variants <- rbind(variants,x)
		}
	}
	return(variants)

} 


is.zero <- function(x){
	print(paste("número:",x))
	if(x == 0){
		return(TRUE)
	}else{
		return(FALSE)
	}
}