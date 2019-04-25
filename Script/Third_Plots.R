library(ggplot2)
library(PerformanceAnalytics)
library(stringr)
library(tidyverse)
library(dplyr)
library(magrittr)

#######################################################################################
# input file 
tips <- read.csv("dataPathways.csv", sep=",", header = TRUE)
length(unique(tips$pathway)) # unique pathway names
sum(tips$count) # total number of records
rm_paths <- c("Metabolic pathways","Biosynthesis of secondary metabolites",
                  "Biosynthesis of antibiotics","Microbial metabolism in diverse environments")
tipsfinal<- filter(tips, pathway!= rm_paths[1], pathway!= rm_paths[2],
                      pathway!= rm_paths[3],pathway!= rm_paths[4]) #remove generic pathway




# Records
tipsmean <- tipsfinal %>% group_by(tissue, health_state, set) %>%
  summarise(Number_pathways = length(set),Records=sum(count),
            Abundance = Records/Number_pathways)
tips_rep <- tipsfinal %>% group_by(tissue, health_state, set, count)%>% 
  summarise(Repetition = length(factor(count))) %>%
  mutate(prop_rep = Repetition/sum(Repetition))
tips_reps <- tips_rep %>% filter(prop_rep>= 0.02) # select only most common



# create a column of proportions
tips_p <- tips %>% group_by(network_code) %>% mutate(prop = 100*count / sum(count)) # original file
tipsfinal_p <- tipsfinal%>% group_by(network_code) %>% mutate(prop = 100*count / sum(count)) # final file

# select the top10 most abundant pathways (by proportion)
tipsfinal_top10 <- tipsfinal_p %>% group_by(tissue, health_state, set) %>% top_n(10,prop)

# unique pathways by network
paths_by_network <- tipsfinal_p %>% group_by(network_code) %>% select (tissue, health_state,set, network_code, count) %>%
  summarise(Tissue =unique(tissue), Health = unique(health_state),set = unique(set), 
            number_of_pathways = length(network_code), Abundance = sum(count))

# unique networks by pathways
networks_by_path <- tipsfinal_p %>% group_by(pathway) %>% summarise(number_of_networks = length(pathway))
network_path_dist <- networks_by_path %>% group_by(number_of_networks) %>% summarise(frec = length(number_of_networks))
network_path_dist <- network_path_dist %>% group_by(number_of_networks, frec)%>% 
  summarise(prop_frec = 100*frec/sum(network_path_dist$frec))
paths_by_network$Ratio <- paths_by_network$Abundance/paths_by_network$number_of_pathways


# Input data for plot 3A: enzyme records in original data
primaryData <- read.csv("HMR_counts_initialData.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
ECs <- primaryData %>% group_by(Tissue, HealthState, Set, Enzyme) %>% summarise(Count = length(Enzyme))
ECs_plot <- ECs%>%group_by(Tissue, HealthState, Set) %>% summarise(Number_ECs = length(Set),
                                                                   Records=sum(Count),
                                                                   Abundance = Records/Number_ECs)

# Input data for plot 6A: pairs of sets and percentage of difference
inputplot <- read.csv("Ftest_ProportionsPlot.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE) 



##############################  PLOTs ##########################################

# 3A enzymes and abundance
ggplot(ECs_plot, aes(y=Number_ECs, x = Set, col=HealthState, size=Abundance))+ 
  geom_point(alpha=0.8)+ facet_grid(.~Tissue)+ scale_size_continuous(range = c(3, 10))+
  ylab("Number of enzyme codes")+ xlab("Tissue ~ Set") +ggtitle("3A")+ theme_bw()

# 3B pathways and abundance 
ggplot(paths_by_network, aes(x = set, y= number_of_pathways, col=Health, size = Abundance)) + 
  geom_point(alpha=0.9)+ theme_bw()+  facet_grid(.~Tissue) + ggtitle("3B") +
  ylab("Number of unique pathways") + xlab("Tissue ~ Set") + expand_limits(y = c(100, 125)) + 
  scale_y_continuous(breaks = c(100,105,110,115,120,125,130)) +
  scale_colour_discrete(name = "Health state") + scale_size(range = c(2, 15))

# 4A top10 pathways ####################
ggplot(data = tipsfinal_top10, aes(y = prop, x = set, fill = pathway)) + 
  geom_bar(stat="identity", width = 0.8, colour = "White") +
  theme(legend.direction ="vertical") + xlab("Tissue ~ Set") + 
  ylab("Pathways proportion (%)") + scale_fill_discrete(guide = guide_legend(ncol =1))+
  scale_y_continuous(breaks=c(5,10, 15, 20, 25, 30,35,40,45,50,55,60))+ 
  facet_grid(health_state~tissue) + theme_bw() + ggtitle("4A")

# 4B boxplot of proportional distribution of abundance
ggplot(tipsfinal_p, aes(y= prop, x= set, col=health_state)) + geom_boxplot()  + theme_bw()+ 
  facet_grid(.~tissue) + ggtitle("4B") + 
  ylab("Proportion (%)") +  xlab("Tissue ~ Set") +  scale_y_log10() +
  scale_colour_discrete(name = "Health state")

# 5A Histogram of number of networks where a pathways is present
ggplot(networks_by_path, aes(x = number_of_networks)) + 
  geom_histogram(stat = "bin") +theme_bw() + ggtitle("5A") +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_x_continuous(breaks = c(2,4,6,8,10,12,14,16))+
  ylab("Number of pathways")+ xlab("Number of networks")

# 6A Comparison of percentage of diference
ggplot(inputplot, aes(x=Set1, y= Set2)) +  # plot percent
  geom_point(aes(size=Difference, fill = Difference), shape = 21) + theme_bw() +
  xlab("Set1") + ylab("Set2") + 
  scale_fill_viridis_c(guide = "legend") +
  scale_size_continuous(range = c(2, 20))

