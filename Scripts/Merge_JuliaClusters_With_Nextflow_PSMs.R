library(dplyr)
setwd("Documents/githubspot/meltmatch/")
#make the file for summarizing the julia membership

full_psms <- read.csv("julia_clustering_input_names.csv")
full_psms <- full_psms$x

membership <- read.csv("analysis/Julia/partition_membership.csv")

membership$scan_name <- gsub("bin[0-9]+__", "", full_psms)
write.csv(membership, "output/julia_clustering_membership_leiden.csv")
rm(full_psms, membership)
gc()

#start from here if its already made
membership <- read.csv("output/julia_clustering_membership_leiden.csv")
membership$scan_name <- gsub("1000_", "", membership$scan_name)
#merge with target psm table
#read from meltome
searched_psms <- read.delim("~/Documents/githubspot/deepmeltome/data/target_psmtable.txt")
#matching info is stored in "SpecID" column
searched_psms$scan_name <- searched_psms$SpecID
searched_psms  <- left_join(searched_psms, membership, by = "scan_name")

#clear
rm(membership)
gc()

#gene IDs per bin
geneID_clusters <- unique(searched_psms[,c("Gene.Name", "cluster")])
write.csv(geneID_clusters, "output/geneID_clusters_membership_leiden.csv")

#plot
geneID_clusters <- geneID_clusters %>%group_by(cluster) %>% summarise(Count = n())



library(ggplot2)
ggplot(geneID_clusters, aes(x = cluster, y = Count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Gene count per cluster", x = "Cluster", y = "Count")




