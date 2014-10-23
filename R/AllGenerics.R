###############################################################################
# Package: pampuCaller
# Script: AllGenerics.R
# 
# Author: smonzon
# date: september 2014
###############################################################################

## getters
setGeneric("get_control",function(object,...) standardGeneric("get_control"))
setGeneric("get_test",function(object,...) standardGeneric("get_test"))
setGeneric("get_samples",function(object,...) standardGeneric("get_samples"))
setGeneric("get_regions",function(object,...) standardGeneric("get_regions"))

## calculations
setGeneric("calling_prep",function(object,...) standardGeneric("calling_prep"))
setGeneric("graph_barerror",function(object,...) standardGeneric("graph_barerror"))
setGeneric("mean_sd",function(object,...) standardGeneric("mean_sd"))
setGeneric("regions_filter",function(object,...) standardGeneric("regions_filter"))

