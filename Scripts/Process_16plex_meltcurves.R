library(dbscan)
library(readr)
library(data.table)


#constants
setwd("~/Documents/GitHub/meltmatch/")
calculated <- TRUE

#calculate from raw first
if(calculated == FALSE){path_to_ms2 <- "data/meltome_all_ms2_quants.txt"
column_meta <- data.frame("temperature" = c(1:8, 8:1),
                          "bins" = c("1", "1", "1", "1", "1", "1", "1", "1",
                                     "2", "2", "2", "2", "2", "2", "2", "2" ))
denominator_col_index <- 1 #set the column non denatured for making fractions

#analysis
#if there are multiple bins, account for 
full_psms <- as.data.frame(read_delim(path_to_ms2))
rownames(full_psms) <- full_psms$spectraID
full_psms$spectraID <- NULL
data_reshape <- as.data.frame(matrix(nrow = 0, ncol = length(unique(column_meta$temperature))))
colnames(data_reshape) <- unique(column_meta$temperature)

for(i in unique(column_meta$bins)) {

  bin_cols <- which(column_meta$bins == i)
  in_loop <- full_psms[, bin_cols]
  

  bin_temps <- column_meta$temperature[bin_cols]
  order_temps <- order(bin_temps)
  

  in_loop <- in_loop[, order_temps]
  colnames(in_loop) <- sort(unique(column_meta$temperature))
  
  rownames(in_loop) <- paste0("bin", i, "__", rownames(in_loop))
  data_reshape <- rbind(data_reshape, in_loop)
  
}

original_channelcols <- colnames(full_psms )
full_psms <- data_reshape
rm(data_reshape, in_loop)
gc()
#Consider that not normalizing would make all rare PSMs extremely far from the source protein due to intensity difference
#We must normalize all as fraction non denatured
denominators <- full_psms[,denominator_col_index ]
full_psms <- full_psms/denominators
saveRDS(full_psms, "adjusted_raw_psms.RDS")}
full_psms <- readRDS("adjusted_raw_psms.RDS")

#remove NA values
full_psms <-full_psms[complete.cases(full_psms),]
fwrite(full_psms, "julia_clustering_input.csv")
write.csv(rownames(full_psms), "julia_clustering_input_names.csv")
gc()
#Consider that you can't match scan by sample/ID
#Each scan and sample first are clustered individually
#Consider two approaches: 1st comprehensive, use tail-down weighted SSN and Leiden graph clustering, second density based if you desire a more robust reproducible melting and removing unusual spectras

#Approach 1
#data too big - use Julia

#library(igraph)
# 
# snn_graph_simple <- function(data, subset,#as a vector
#                                k_neighbors = 5, temp_vector = NULL, id_col = "sampleID", id_weight_min = 0.5, temp_weight_min = 0) {
# 
#   rownames_dat <- rownames(data)
#   subset_data <- data[rownames_dat == subset, ]
#   
#   #Euclidean distance matrix
#   #Function to calculate distance matrix for a chunk
#   calculate_chunk_distance <- function(data_chunk) {
#     dist_matrix_chunk <- as.matrix(dist(data_chunk, method = "euclidean"))
#     return(dist_matrix_chunk)
#   }
#   
#   #Number of rows in each chunk
#   chunk_size <- 5000  
#   
#   #Number of chunks
#   num_chunks <- ceiling(nrow(subset_data) / chunk_size)
#   
#   #Initialize a matrix to store the complete distance matrix
#   dist_matrix_complete <- big.matrix(nrow = nrow(subset_data), ncol = nrow(subset_data), 
#                                                 type = 'integer', init = 5)
#   
#   #Loop over chunks
#   for (i in 1:num_chunks) {
#     for (j in i:num_chunks) {
#       # Define row indices for chunks i and j
#       start_row_i <- ((i - 1) * chunk_size) + 1
#       end_row_i <- min(i * chunk_size, nrow(subset_data))
#       start_row_j <- ((j - 1) * chunk_size) + 1
#       end_row_j <- min(j * chunk_size, nrow(subset_data))
#       
#       # Extract chunks
#       data_chunk_i <- subset_data[start_row_i:end_row_i, ]
#       data_chunk_j <- subset_data[start_row_j:end_row_j, ]
#       
#       # Calculate distance matrix between chunks
#       dist_ij <- as.matrix(dist(rbind(data_chunk_i, data_chunk_j), method = "euclidean"))
#       
#       # Fill the appropriate section of the complete matrix
#       dist_matrix_complete[start_row_i:end_row_i, start_row_j:end_row_j] <- dist_ij[1:(end_row_i - start_row_i + 1), (end_row_i - start_row_i + 1 + 1):ncol(dist_ij)]
#       if (i != j) {
#         dist_matrix_complete[start_row_j:end_row_j, start_row_i:end_row_i] <- t(dist_ij[1:(end_row_i - start_row_i + 1), (end_row_i - start_row_i + 1 + 1):ncol(dist_ij)])
#       }
#     }
#   }
#   
#   #Temperature-based weight matrix
#   # if(is.null(temp_vector)){temp_vector <- c(1:length(colnames(subset_data)))}
#   # temp_weight_matrix <- outer(temp_vector, temp_vector, FUN = function(x, y) ifelse(x < y, sqrt((y - x) / (max(temp_vector) - min(temp_vector))), temp_weight_min))
#   # 
#   
#   # temp_weight_matrix <- outer(temp_vector, temp_vector, FUN = function(x, y) {
#   #   weight = sqrt((y - x) / (max(temp_vector) - min(temp_vector)))
#   #   scaled_weight = (weight * (1 - temp_weight_min)) + temp_weight_min
#   #   return(ifelse(x < y, scaled_weight, temp_weight_min))
#   # })
#   
#   #Sample ID similarity weight matrix
#   sample_id_vector <-  gsub("bin2__", "", gsub("bin1__", "", rownames(subset_data)))
#   id_weight_matrix <- outer(sample_id_vector, sample_id_vector, FUN = function(x, y) ifelse(x == y, 1, id_weight_min))
#   
#   #Combine weight matrices
#   combined_weight_matrix <- #temp_weight_matrix * 
#     id_weight_matrix
#   
#   #Finding k-nearest neighbors for each sample and apply weights
#   knn_graph <- matrix(0, nrow = nrow(dist_matrix_complete), ncol = nrow(dist_matrix_complete))
#   for (i in 1:nrow(dist_matrix_complete)) {
#     neighbors <- sort(dist_matrix_complete[i,], decreasing = FALSE)[1:k_neighbors]
#     knn_graph[i, names(neighbors)] <- combined_weight_matrix[i, names(neighbors)]
#   }
#   
#   # Constructing the weighted SNN graph
#   snn_graph <- graph_from_adjacency_matrix(knn_graph, mode = "undirected", weighted = TRUE, diag = FALSE)
#   
#   return(snn_graph)
# }
# 
# 
# 
# snn_graph_taildown <- function(data, subset, # as a vector
#                                k_neighbors = 5, temp_vector = NULL, 
#                                id_weight_min = 0.5, temp_weight_min = 0) {
#   
#   # Selecting the subset of data
#   subset_data <- data[rownames(data) %in% subset, ]
#   
#   # Temperature-based weight matrix
#   if(is.null(temp_vector)) {
#     temp_vector <- as.numeric(colnames(subset_data)) # Convert column names to numeric
#   }
#   temp_range <- max(temp_vector) - min(temp_vector)
#   temp_weight_matrix <-  outer(temp_vector, temp_vector, FUN = function(x, y) {
#     diff = y - x
#     scaled_diff = diff / temp_range
#     weight = ifelse(scaled_diff > 0, sqrt(scaled_diff), 0) # Use ifelse for vectorized operation
#     scaled_weight = (weight * (1 - temp_weight_min)) + temp_weight_min
#     return(scaled_weight)
#   })
#   
#   # Sample ID similarity weight matrix
#   sample_id_vector <- rownames(subset_data)
#   id_weight_matrix <- matrix(id_weight_min, nrow = length(sample_id_vector), ncol = length(sample_id_vector))
#   diag(id_weight_matrix) <- 1
#   
#   # Check if dimensions match and combine weight matrices
#   if(all(dim(temp_weight_matrix) == dim(id_weight_matrix))){
#     combined_weight_matrix <- temp_weight_matrix * id_weight_matrix
#   } else {
#     stop("The dimensions of temperature and ID weight matrices do not match.")
#   }
#   
#   # Finding k-nearest neighbors for each sample and apply weights
#   knn_graph <- matrix(0, nrow = nrow(dist_matrix), ncol = nrow(dist_matrix))
#   for (i in 1:nrow(dist_matrix)) {
#     neighbors <- sort(dist_matrix[i, ], decreasing = FALSE)[1:k_neighbors]
#     neighbor_names <- names(neighbors)
#     if(all(neighbor_names %in% colnames(combined_weight_matrix))){
#       knn_graph[i, neighbor_names] <- combined_weight_matrix[i, neighbor_names]
#     } else {
#       stop("Neighbor names do not match with weight matrix column names.")
#     }
#   }
#   
#   # Constructing the weighted SNN graph
#   snn_graph <- graph_from_adjacency_matrix(knn_graph, mode = "undirected", weighted = TRUE, diag = FALSE)
#   
#   return(snn_graph)
# }
# 
# 
# psms_graph <- snn_graph_simple(data = full_psms, subset = rownames(full_psms), id_weight_min = 1)




#Approach 2
#WAY too slow


#full_psms <- as.big.matrix(full_psms)
# dbscan_clusters <- dbscan::hdbscan(full_psms, gen_hdbscan_tree = F, minPts = 5)
# saveRDS(dbscan_clusters, "dbscan_clusters.RDS")
