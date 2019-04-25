# Packages
library(KEGGREST)
library(stringr)
library(tidyverse)
library(dplyr)

###############################################################################################
################################### Pathways function #########################################
Pathways <- function (file){
  setwd("Data/")
  filename<- paste0('EC_',file,'.txt')
  inputfile <- read.delim(filename, header = FALSE, 
                          sep = "\t", stringsAsFactors = FALSE)
  enzymes <- paste0("ec:",inputfile[,3]) # select only enzymes and write in EC number format 
  data<- vector("list", length(enzymes)) # create an empty vector
  
  # extract all information from KEGG using KeggGet() function from KEGGREST package
  for(i in 1:length(enzymes)){
    data[i] <- tryCatch(keggGet(enzymes[i]), error=function(e) NULL)
  }
  # from "data" select information of pathways  
  l_pathways <- lapply (data, function(m) m["PATHWAY"]) # extract list of pathways 
  df_pathways <- data.frame(pathway_name= unlist(l_pathways)) # tranform in data frame
  c_pathways <- df_pathways %>% group_by(pathway_name) %>%
    summarise(count = length(pathway_name))# count pathway ocurrences
  
  # output files  
  l_outputname <- paste0('LP_EC_',file,'.csv') # LP: list of all pathways
  c_outputname <- paste0('CP_EC_',file,'.csv') # CP: Count of ocurrence for each pathways
  write.csv(df_pathways, file = l_outputname)# list of all pathways
  write.csv(c_pathways, file = c_outputname)# Pathways with nnumber of ocurrences
}

################################### Run Pathways function ##################################### 
# Run "Pathways" for multiple files
### Tissue: Breast, Lung, (Kidney or Renal), Urothelial
### Condition: Cancer, Healthy
### Set: MDS, nonMDS; where MDs: Minimum Dominant Set

names <- data.frame(files= c("BreastCancer_MDS","BreastCancer_nonMDS",
           "BreastHealthy_MDS","BreastHealthy_nonMDS",
           "LungCancer_MDS","LungCancer_nonMDS",
           "LungHealthy_MDS","LungHealthy_nonMDS",
           "RenalCancer_MDS","RenalCancer_nonMDS", 
           "KidneyHealthy_MDS", "KidneyHealthy_nonMDS",
           "UrothelialCancer_MDS","UrothelialCancer_nonMDS",
           "UrothelialHealthy_MDS","UrothelialHealthy_nonMDS"), 
           code= c("BCM","BCnM", "BHM", "BHnM",
                   "LCM","LCnM", "LHM", "LHnM",
                   "RCM","RCnM", "KHM", "KHnM",
                   "UCM","UCnM", "UHM", "UHnM"))
trialnames <- data.frame(files= c("trial1", "trial2", "trial3", "trial4"),
                code = c("T1", "T2", "T3", "T4"))

for(name in names$files){ # change trialnames$files by names$files
  Pathways(name)
}

# Run "Pathways" for single file
Pathways("trial1") # enter name, example "trial1", "BreastCancer_MDS"


###############################################################################################
################################### Statistics function #######################################
statistical_test <- function (count1, count2, testing){
  # Files to be used
  name1<- paste0('Pathway_count/CP_EC_',count1,'.csv')
  input1_temp <- read.csv(name1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  name2<- paste0('Pathway_count/CP_EC_',count2,'.csv')
  input2_temp <- read.csv(name2, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  # remove generic pathways
  rm_paths <- c("Metabolic pathways","Biosynthesis of secondary metabolites",
                "Biosynthesis of antibiotics","Microbial metabolism in diverse environments")
  input1 <- filter(input1_temp, pathway_name!= rm_paths[1], pathway_name!= rm_paths[2],
                   pathway_name!= rm_paths[3],pathway_name!= rm_paths[4]) #remove generic pathway
  input2 <- filter(input2_temp, pathway_name!= rm_paths[1], pathway_name!= rm_paths[2],
                   pathway_name!= rm_paths[3],pathway_name!= rm_paths[4]) #remove generic pathway
  
  # Input values for test
  if(testing == "pt"){
    col1name<- paste0('ptest_',count1,'_', count2) # P: prop.test()
    } else {
      col1name<- paste0('ftest_',count1,'_', count2) # F: fisher.test()
    }
  # Variables to test
  stat_test <- data.frame(file_name = col1name, pathway_name = unique(input1$pathway_name))
  for (nam in stat_test$pathway_name){
    a1 <- input1$count[input1$pathway_name == nam] # first x (simple) value
    b1 <- input2$count[input2$pathway_name == nam] # second x (simple) value
    at <- sum(input1$count) # first n (total) value
    bt <- sum(input2$count) # second n (total) value
    if(is_empty(a1) | is_empty(b1)){
      next
      } else if(testing == "pt"){ 
      stat_test$p_proptest[stat_test$pathway_name==nam] <- prop.test(c(a1, b1), c(at, bt))$p.value
    
      } else {
      stat_test$p_fishertest[stat_test$pathway_name==nam] <- fisher.test(matrix(c(a1, at, b1, bt), 
                                                                           nrow = 2))$p.value
      }
  }
  if(testing == "pt"){
    statname<- paste0('Statistics_PropTest/ptest_',count1,'_', count2, '.csv') # P: prop.test()
    #p_outdata <- stat_test %>% select(file_name, pathway_name,p_proptest)
    write.csv(stat_test, file = statname)
  } else {
    statname<- paste0('Statistics_FisherTest/ftest_',count1,'_', count2, '.csv') # F: fisher.test()
    #f_outdata <- stat_test %>% select(file_name, pathway_name,p_fishertest)
    write.csv(stat_test , file = statname)
  }
}

################################### Run Statistics function ##################################### 
# prop.test() is suggested to compare two proportions using count data (integer values)
# fisher.test() is suggested for two nominal variables evaluated for independence.
# fisher.test() accepts integer values, it is arranged in a contingency table. 
# Files: "file1", "file2"
# Statistical test: "pt" for prop.test(), "ft" for fisher.test()

# Run Statistics function for all files
# Combination of pairs for ptest() or ftest()
comb_names <- as.vector(names$files)
combs <- data.frame(combn(comb_names, 2)) # 
comb1 <- t(combs[1,])
comb2 <- t(combs[2,])

for (i in 1:nrow(comb1)){
  statistical_test(comb1[i], comb2[i], "ft")# pt: prop.test(), ft: fisher.test()
}

#########  Run Statistics function for two single files
# change folder name within function to 
# 'Statistics_FisherTest/ftest_' or 'Statistics_PropTest/ptest_'
statistical_test("trial1", "trial2", "pt") 

#################################################################################################
#################### Function join files and extract p-values <= 0.05 ###########################
# join all files ptest or ftest files
library(magrittr)

join_function <- function (testing){
  if (testing == "pt"){
    setwd("Statistics_PropTest/") # change to "../Statistics_FisherTest/"
    ptest_files <- list.files()
    combined_ptest <- do.call(rbind, lapply(ptest_files, function(x) read.csv(x, header = TRUE, stringsAsFactors = FALSE)))
    psignificant_list <- combined_ptest %>% group_by(file_name)%>%subset(p_proptest <= 0.05)
    psignificant_count <- psignificant_list %>% group_by(file_name)%>%summarise(n_pathways = length(file_name))
    write.csv(combined_ptest, file = "../Ptest_joined.csv") # change to "../Ptest_joined.csv"
    write.csv(psignificant_list , file = "../Ptest_significantList.csv") #
    write.csv(psignificant_count , file = "../Ptest_significantCount.csv") #
  } else if (testing == "ft"){
    setwd("Statistics_FisherTest/") # change to "../Statistics_FisherTest/"
    ftest_files <- list.files()
    combined_ftest <- do.call(rbind, lapply(ftest_files, function(x) read.csv(x, header = TRUE, stringsAsFactors = FALSE)))
    combined_ftest_count <- combined_ftest %>% group_by(file_name)%>% summarise(n_pathways = length(file_name))
    fsignificant_list <- combined_ftest %>% filter(p_fishertest <= 0.05)
    fsignificant_count <- fsignificant_list %>% group_by(file_name)%>%summarise(n_pathways = length(file_name))
    write.csv(combined_ftest, file = "../Ftest_joined.csv") # change to "../Ptest_joined.csv"
    write.csv(combined_ftest_count, file = "../Ftest_joinedCount.csv") # change to "../Ptest_joined.csv"
    write.csv(fsignificant_list , file = "../Ftest_significantList.csv") #
    write.csv(fsignificant_count , file = "../Ftest_significantCount.csv") #
    }
  setwd("../") # return to main file
}
############################# Run join ans p-value <= 0.05 function ################################# 
join_function ("ft") # use "ft" to join all fisher test files

