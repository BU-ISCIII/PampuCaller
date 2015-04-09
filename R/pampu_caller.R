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
#CHR=as.character(),POS=as.numeric(),REF=as.character(),A=as.numeric(),C=as.numeric(),T=as.numeric(),
		#G=as.numeric(),a=as.numeric(),c=as.numeric(),t=as.numeric(),g=as.numeric/
pampu_caller_test_1 <- function(x,control,min_freq,times){
	variants <- data.frame(CHR=as.character(),POS=as.numeric(),REF=as.character(),A=as.numeric(),C=as.numeric(),T=as.numeric(),
		G=as.numeric(),a=as.numeric(),c=as.numeric(),t=as.numeric(),g=as.numeric(),sample=as.character(),depth=as.numeric(),
		total_A=as.numeric(),total_C=as.numeric(),total_T=as.numeric(),total_G=as.numeric(),forward_reads=as.numeric(),
		reverse_reads=as.numeric(),per_forward_reads=as.numeric(),per_reverse_reads=as.numeric(),per_A=as.numeric(),per_C=as.numeric(),
		per_T=as.numeric(),per_G=as.numeric(),per_a=as.numeric(),per_c=as.numeric(),per_t=as.numeric(),per_g=as.numeric(),
		per_total_A=as.numeric(),per_total_C=as.numeric(),per_total_T=as.numeric(),per_total_G=as.numeric(),ALT=as.character(),stringsAsFactors=FALSE)

	print(class(x))
	print(x)
	seek_control_pos <- control[control["POS"] == as.numeric(x["POS"]),]
	print(seek_control_pos)

	#print(paste("POS:",x["POS"],"testing"))
	if(nrow(seek_control_pos["mean_per_total_A"]) == 0 || nrow(seek_control_pos["mean_per_total_C"]) == 0 || nrow(seek_control_pos["mean_per_total_G"])==0 || nrow(seek_control_pos["mean_per_total_T"])==0){ 
		print(paste("POS:",x["POS"],"not found in control."))
		return(variants)
	}

	## Test each nt independently 
	if(is.zero(as.numeric(seek_control_pos["mean_per_total_A"]))){
		if(x["per_total_A"] > as.numeric(times/as.numeric(x["depth"]))){
			x["ALT"] <- "A"
			variants <- rbind(variants,x)
			print("llega")
		}
	}else{
		if(x["per_total_A"] > times*as.numeric(seek_control_pos["mean_per_total_A"])){
			x["ALT"] <- "A"
			variants <- rbind(variants,x)
			print("llega")
		}
	}
	if(is.zero(as.numeric(seek_control_pos["mean_per_total_C"]))){
		if(x["per_total_C"] > as.numeric(times/as.numeric(x["depth"]))){
			x["ALT"] <- "C"
			variants <- rbind(variants,x)
			print("llega")
		}
	}else{
		if(x["per_total_C"] > times*as.numeric(seek_control_pos["mean_per_total_C"])){
			x["ALT"] <- "C"
			variants <- rbind(variants,x)
			print("llega")
		}
	}
	if(is.zero(as.numeric(seek_control_pos["mean_per_total_T"]))){
		if(x["per_total_T"] > as.numeric(times/as.numeric(x["depth"]))){
			x["ALT"] <- "T"
			print(x)
			variants <- as.data.frame(rbind(variants,as.data.frame(x)))
			print(variants)
			print("llega")
		}
	}else{
		if(x["per_total_T"] > times*as.numeric(seek_control_pos["mean_per_total_T"])){
			x["ALT"] <- "T"
			print(x)
			variants <- as.data.frame(rbind(variants,as.data.frame(x)))
			print(variants)
			print("llega")
		}
	}
	if(is.zero(as.numeric(seek_control_pos["mean_per_total_G"]))){
		if(x["per_total_G"] > as.numeric(times/as.numeric(x["depth"]))){
			x["ALT"] <- "G"
			print(x)
			variants <- as.data.frame(rbind(variants,x))
			print(variants)
			print("llega")
		}
	}else{
		if(x["per_total_G"] > times*as.numeric(seek_control_pos["mean_per_total_G"])){
			x["ALT"] <- "G"
			print(x)
			variants <- as.data.frame(rbind(variants,as.data.frame(x)))
			print(variants)
			print("llega")
		}
	}
	
	return(variants)

} 


is.zero <- function(x){
	if(x == 0){
		return(TRUE)
	}else{
		return(FALSE)
	}
}