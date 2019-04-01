# Packages
library(KEGGREST)
library(stringr)
library(tidyverse)
library(dplyr)

###############################################################################################
################################### Pathways function #########################################
Pathways <- function (file){
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
# Run "Pathways" for single file
Pathways("BreastHealthy_nonMDS") # enter name, example "trial1", "BreastCancer_MDS"

# Run "Pathways" for multiple files
### Tissue: Breast, Lung, Kidney, Renal, Urothelial
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

for(name in names$files){ # change trialnames by names$files
  Pathways(name)
}

###############################################################################################
################################### Statistics function #######################################
statistical_test <- function (count1, count2, testing){
  # Files to be used
  name1<- paste0('CP_EC_',count1,'.csv')
  input1 <- read.csv(name1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  name2<- paste0('CP_EC_',count2,'.csv')
  input2 <- read.csv(name2, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  
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
    statname<- paste0('ptest_',count1,'_', count2, '.csv') # P: prop.test()
    #p_outdata <- stat_test %>% select(file_name, pathway_name,p_proptest)
    write.csv(stat_test, file = statname)
  } else {
    statname<- paste0('ftest_',count1,'_', count2, '.csv') # F: fisher.test()
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

# Run Statistics function for two single files
statistical_test("trial1", "trial2", "pt")

# Run Statistics function for all files
test_option <- c("pt", "ft") # pt: prop.test(), ft: fisher.test()
for(name1 in trialnames){ # change trialnames by names
  for(name2 in trialnames){
    for(tst in test_option){
      statistical_test(name1, name2, tst)
    }
  }
}


#################################################################################################
############################## Significant p-value function #####################################

# join all files ptest or ftest files
ptest_files <- list.files(pattern = "ptest_*")
ptest_allfiles = do.call(rbind, lapply(ptest_files, function(x) read.csv(x, stringsAsFactors = FALSE)))


# significant p-values
library(magrittr)
ptest_significant <- ptest_allfiles %>% group_by(file_name)%>%
  subset(p_proptest <= 0.05) %>% summarise(n_pathways = length(file_name))

pdf("Important_pathways.pdf")
ggplot(ptest_significant, aes(x=file_name, y= n_pathways)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Pathways significantly important") +
  xlab("Compared data") +
  ylab("Pathways (n)")
dev.off()















##############################################################################################
######## Fisher test and Proportional test comparison (linear model and plot) ################

statname <- paste0('stat_',filenames[1],'_', filenames[2], '.csv')
statfile <- read.csv(statname, header = TRUE, sep = ",", stringsAsFactors = FALSE)
linearMod <- lm(fishertest ~ proptest, data=statfile)
summary(linearMod)

####################### Plot Linear Regression Function ######################################
# This function is to validate the that prop.test and fisher.test p-values are similar
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() + geom_jitter() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

####################### Run Linear Regression Function ########################################
plotname <- paste0(filenames[1],'_', filenames[2], '.jpg')
jpeg(plotname, width = 960, height = 960)
ggplotRegression(linearMod) # run only this line to visualise plot
dev.off()

####################### Logistic Regression Function ##########################################
logitMod <- glm(fishertest ~ proptest, data=statfile, family=binomial(link="logit"))
