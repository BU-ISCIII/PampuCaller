cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

brewer_qualitative <- c("#0000ff","#ff0000","#483215","#008900","#7244c4","#e65a11","#000000","#e6e528","#ff00ee","#6e0000","#00c7dd","#d4b455","#8f008d","#736b00","#7d8cbf")

################################################
### PRUEBA Octubre 2014
############ ####################################

library(pampuCaller)
targets <- read_targets()
regions <- read_regions("amplicones.bed")
position_set <- create_position_set(targets,regions) 
compare_control_set <- calling_prep(position_set,depth=100)

dir.create("graphs_R")
apply(as.matrix(get_regions(compare_control_set)),1,draw_graph,"graphs_R",set=compare_control_set)


draw_graph <- function(reg,dir,set){
		print(paste("Graph",reg[4]))
		
		dir <- paste(dir,"/",reg[4],sep="")

		dir.create(dir)
		test <- get_test(compare_control_set)
		pos_reg <- pos_filter(test,reg,depth=0)
		depth_stats <- ddply(pos_reg,.(sample),summarize,mean_depth=mean(depth))
		write.table(depth_stats,file=paste(dir,"/",reg[4],"_depth_stats.txt",sep=""),row.names=FALSE,sep="\t")
		g <- graph_barerror(set,start_end=reg,zoom=c(0,0.015))
		pdf(paste(dir,"/",reg[4],"_zoom_0.015",sep=""),width=80)
			print(g[[1]])
			print(g[[2]])
		dev.off()

		g <- graph_barerror(set,start_end=reg,zoom=c(0,0.008))
		pdf(paste(dir,"/",reg[4],"_zoom_0.008",sep=""),width=80)
			print(g[[1]])
			print(g[[2]])
		dev.off()

		g <- graph_barerror(set,start_end=reg,zoom=c(0,0.5))
		pdf(paste(dir,"/",reg[4],"_zoom_0.5",sep=""),width=80)
			print(g[[1]])
			print(g[[2]])
		dev.off()

		g <- graph_barerror(set,start_end=reg,zoom=c(0,0.25))
		pdf(paste(dir,"/",reg[4],"zoom_0.25",sep=""),width=80)
			print(g[[1]])
			print(g[[2]])
		dev.off()

		g <- graph_barerror(set,start_end=reg,zoom=c(0,1))
		pdf(paste(dir,"/",reg[4],"zoom_1.0",sep=""),width=80)
			print(g[[1]])
			print(g[[2]])
		dev.off()
		g <- graph_barerror(set,start_end=reg,zoom=c(0,0.75))
		pdf(paste(dir,"/",reg[4],"zoom_1.0",sep=""),width=80)
			print(g[[1]])
			print(g[[2]])
		dev.off()
}



## Test
position_set_filt <- regions_filter(position_set,depth=1000)
pfc <- get_control(position_set_filt)
a <- mean_sd(position_set_filt)

start_end <- get_regions(compare_control_set)[1,]
depth <- 0
zoom <- c(0,0.008)
samples <- compare_control_set@samples$Sample
var <- get_test(compare_control_set)
var_c <- get_control(compare_control_set)

l <- graph_barerror(compare_control_set,start_end=start_end,zoom=c(0,0.008))

## Obtain data from specific mutations.
mut <- targets
mut$merge <- paste(mut$Sample,"_",mut$pos,sep="")
test <- get_test(position_set)
test$merge <- paste(test$sample,"_",test$POS,sep="")
test_validation <- test[test$merge %in% mut$merge,]
merge_mut_validation <- merge(mut,test_validation,by="merge")
write.table(merge_mut_validation,file="mutations_validation_deepseq.txt",row.names=F,sep="\t")


test <- get_test(compare_control_set)
control <- get_control(compare_control_set)
x <- test[test$POS == 156713 & test$sample=="3543_rep1",]
p <- apply(x,1,pampu_caller_test_1,control=control,times=10)
