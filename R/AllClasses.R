###############################################################################
# Package: pampuCaller
# Script: AllClasses.R
# 
# Author: smonzon
# date: september 2014
###############################################################################


## controlSet Class.

setClass("PositionSet",
		representation(control = "data.frame",
					   test="data.frame",
					   regions = "data.frame",
					   polymorphisms = "data.frame",
					   samples = "data.frame"
						)
		)

setClass("ControlCompareSet",
		representation(meanControl = "data.frame",
					   test="data.frame",
					   regions = "data.frame",
					   polymorphisms = "data.frame",
					   samples = "data.frame"
						)
		)

setClass("VariationSet",
		representation(variants="data.frame",
					   annotation = "data.frame"
						)
		)