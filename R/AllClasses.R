###############################################################################
# Package: pampuCaller
# Script: AllClasses.R
# 
# Author: smonzon
###############################################################################


## controlSet Class.

setClass("PositionSet",
		representation(variation = "data.frame",
						genes = "data.frame",
						polymorphisms = "data.frame",
						samples = "data.frame"
						)
		)