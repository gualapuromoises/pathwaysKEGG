library(tidyverse)
library(magrittr)
library(ggplot2)
library(tidyr)

# data file

tips <- read.csv("dataPathways.csv", sep=",", header = TRUE)
length(unique(tips$pathway)) # unique pathway names
sum(tips$count) # total number of records
rm_paths <- c("Metabolic pathways","Biosynthesis of secondary metabolites",
              "Biosynthesis of antibiotics","Microbial metabolism in diverse environments")
dat<- filter(tips, pathway!= rm_paths[1], pathway!= rm_paths[2],
                   pathway!= rm_paths[3],pathway!= rm_paths[4]) #remove generic pathway

# FISHER EXACT TEST
# Suggested for two nominal variables evaluated for independence.
# This test accepts integer values, because it works with count data arranged in a contingency table. 
fisher_test <- data.frame(path = unique(dat$pathway))
for (nam in fisher_test$path){
  a1 <- dat$count[dat$pathway == nam & dat$health_state=="Cancer"& dat$tissue=="Urothelial"]
  b1 <- dat$count[dat$pathway == nam & dat$health_state=="Healthy"& dat$tissue=="Urothelial"]
  c1 <- sum(dat$count[dat$health_state=="Cancer"& dat$tissue=="Urothelial"])
  d1 <- sum(dat$count[dat$health_state=="Cancer"& dat$tissue=="Urothelial"])
  if(is_empty(a1) | is_empty(b1) | is_empty(c1) | is_empty(d1)){
    next
  } else{
    fisher_test$pvalue[fisher_test$path==nam] <- fisher.test(matrix(c(a1, c1, b1, d1), nrow = 2))$p.value
  }
}
write.csv(fisher_test, file="UrothelialF_CH.csv")


# FISHER EXACT TEST (trial with proportions)
# Double values are ROUNDED to integer values. Does not work properly  
fisherp_test <- data.frame(pathway = unique(dat$pathway))
for (nam in fisher_test$pathway){
  a <- dat$count[dat$pathway == nam & dat$network_code=="UCM"]
  b <- dat$count[dat$pathway == nam & dat$network_code=="UcnM"]
  c <- dat$count[dat$pathway == nam & dat$network_code=="UHM"]
  d <- dat$count[dat$pathway == nam & dat$network_code=="UhnM"]
  if(is_empty(a) | is_empty(b) | is_empty(c) | is_empty(d)){
    next
  } else{
    a1 <- a/sum(dat$Count[dat$network_code=="UCM"])
    b1 <- b/sum(dat$Count[dat$network_code=="UcnM"])
    c1 <- c/sum(dat$Count[dat$network_code=="UHM"])
    d1 <- d/sum(dat$Count[dat$network_code=="UhnM"])
    fisherp_test$Urothelialp[fisherp_test$path_name==nam] <- fisher.test(matrix(c(a1, c1, b1, d1), nrow = 2))$p.value
  }
}
write.csv(fisherp_test, file="fisher_proportions.csv")

# plots for overview
ggplot(dat, aes(x=network_code, y=count, color=tissue)) + geom_violin(alpha=0.2)+ geom_point(alpha=0.7)+
  scale_y_log10()

ggplot(dat, aes(count, color=health_state)) + facet_grid(.~tissue)+ geom_density(kernel = "gaussian")+ scale_x_log10()
