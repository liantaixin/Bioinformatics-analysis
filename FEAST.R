
rm(list = ls())

library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("reshape2")
library("ggplot2")
library("cowplot")
library("Rcpp")
library("RcppArmadillo")
library("FEAST")


source("src.R")
metadata_file = "groupHem.csv" #group
count_matrix = "ASVs-CON.csv" #otu table
metadata <- read.csv(file = metadata_file, header=T, sep = ",", row.names = 1)
otus <- read.csv(file = count_matrix, header=T, sep = ",", row.names = 1)
otus <- t(as.matrix(otus))
EM_iterations = 1000 #default 1000
different_sources_flag = 0
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
if(length(common.sample.ids) <= 1) {  # 如果common.sample.ids的长度小于等�?1  
  message <- paste(sprintf('Error: there are %d sample ids in common '),     
  'between the metadata file and data table') 
# message 是一个字符串，表示错误信息，显示common.sample.ids的长�? 
stop(message) 
# 停止运行，并显示message
}
if(different_sources_flag == 0){
  # 如果different_sources_flag等于0
  metadata$id[metadata$SourceSink == 'Source'] = NA
  # 将metadata中SourceSink列等�?'Source'的id列赋值为NA（缺失值）
  metadata$id[metadata$SourceSink == 'Sink'] = c(1:length(which(metadata$SourceSink == 'Sink')))
  # 将metadata中SourceSink列等�?'Sink'的id列赋值为1到SourceSink列等�?'Sink'的个数}
}
  envs <- metadata$Env
  Ids <- na.omit(unique(metadata$id))
  Proportions_est <- list()  
  for(it in 1:length(Ids)){
    if(different_sources_flag == 1){
      train.ix <- which(metadata$SourceSink=='Source' & metadata$id == Ids[it])
      test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])
    } else {
      train.ix <- which(metadata$SourceSink=='Source')
      test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])
    }
    num_sources <- length(train.ix)
    COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),])) 
    str(COVERAGE)
    sources <- as.data.frame(as.matrix(rarefy(as.data.frame(otus[train.ix,]), COVERAGE)))
    sinks <- as.data.frame(as.matrix(rarefy(as.data.frame(t(as.matrix(otus[test.ix,]))), COVERAGE)))
    print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
    print(paste("Seq depth in the sources and sink samples = ",COVERAGE))
    print(paste("The sink is:", envs[test.ix]))
    FEAST_output<-FEAST(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
    Proportions_est[[it]] <- FEAST_output$data_prop[,1]
    names(Proportions_est[[it]]) <- c(as.character(envs[train.ix]), "unknown")
    if(length(Proportions_est[[it]]) < num_sources +1){
      tmp = Proportions_est[[it]]
      Proportions_est[[it]][num_sources] = NA
      Proportions_est[[it]][num_sources+1] = tmp[num_sources]
    }
    print("Source mixing proportions")
    print(Proportions_est[[it]])
  }
  write.csv(Proportions_est,file = "res2.csv",quote=F)
  
  