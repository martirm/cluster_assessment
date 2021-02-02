library(igraph)
library(igraphdata)
set.seed(1)

### CLUSTERING ALGORITHMS
cluster_leading_eigen_tryCatch <- function(g){
    c <- tryCatch({
        #do.call(f, g)
        cluster_leading_eigen(g)
    }, 
    error=function(c){
        print("error")
        message(c)
        return(make_clusters(g, replicate(gorder(g), 1)))
    }
    # ,
    # warning=function(c){
    #     print("warning")
    #     message(c)
    #     return(make_clusters(g, replicate(gorder(g), 1)))
    # }
    )
    return(c)
}

non_connected_spinglass <- function(g){
    # Splits the graph into connected components, clusters each one with the 
    # spinglass algorithm, and then merges them again (reindexing the 
    # communities appropiately).
    component_list <- decompose.graph(g)
    memb <- function(g_component) {
        if (gorder(g_component)>1){
            return(membership(cluster_spinglass(g_component)))
        }
        else{
            #solution to avoid crashes when g_component is an isolated vertex
            m <- as.vector(1)
            names(m) <- names(V(g_component))
            return(m)
        }
        
    }
    membership_list <- lapply(component_list, memb)
    count_comms <- function(m) length(unique(m))
    n_comms <- lapply(membership_list, count_comms)
    offset <- n_comms %>% cumsum %>% c(0, .) %>% head(-1) %>% as.list
    combined_membership <- mapply("+", membership_list, offset) %>% unlist
    make_clusters(g, membership = combined_membership)
}

clust_alg_list <- c(cluster_louvain, cluster_leading_eigen_tryCatch, cluster_label_prop, cluster_walktrap, non_connected_spinglass)
names(clust_alg_list) <- c("Louvain", "leading eigen", "label prop", "Walktrap", "spinglass")


### GRAPHS
#news graph
id_dict <- read.csv(file="data/id_dict.csv", sep="|", stringsAsFactors=FALSE, encoding = "UTF-8")
edge_list <- read.csv(file="data/news_graph.csv")
g_news <- graph_from_data_frame(edge_list, directed=FALSE)
for (i in V(g_news)){
    V(g_news)[i]$name <- id_dict[V(g_news)[i],2]
    #set_vertex_attr(g_news, name="name", index=i, value=id_dict[i,2])
}
g_news_small <- induced_subgraph(g_news,degree(g_news)>100)

#Enron
data(enron)
weighted_enron <- graph.adjacency(get.adjacency(as.undirected(enron, mode="each")),weighted=TRUE, mode="undirected")
#E(weighted_enron)$weight <- 1
#weighted_enron <- simplify(weighted_enron, edge.attr.comb=list(weight="sum", "ignore"))
V(weighted_enron)$name <- V(enron)$Name
#weighted_enron <- weighted_enron %>% delete_vertex_attr("Name") %>%
#    delete_vertex_attr("Note") %>% delete_vertex_attr("Email")


#Opsahl social network
#Facebook-like social network in https://toreopsahl.com/datasets/#online_forum_network (first dataset)
#vertices are users, edges are the character count of their interactions.
social_network <- read.csv(file="data/Opsahl_social_network_data.txt", sep=" ", stringsAsFactors=FALSE)
social_network <- social_network[social_network[,2]!=social_network[,3],]
g_social_network <- graph_from_data_frame(social_network[,2:4], directed=FALSE) #network is actually directed, but listed algorithms are not compatible with that (look into this!)
#g_social_network <- decompose.graph(g_social_network)[[1]] #this removes 3 isolated pairs of vertices 
# (which allows the use of algorithms with implementations incompatible with disconnected graphs)

#CondMat graph
#CondMat <- read.csv(file="data/CA-CondMat.txt", sep=" ", stringsAsFactors=FALSE)

#Zachary's karate club
library(igraphdata)
data(karate)
karate_gt_clustering <- c(1,1,1,1,1,1,1,1,1,2,1,1,1,1,2,2,1,1,
                          2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)
#corresponds to the final club of each member after the split

#Forex data
# source("main_new.R")
# g_forex <- graph_returns(fx_returns())
load("data/g_forex.RData")

#Stochastic block model with 4 clusters of 25 nodes
#matrix of probabilites of edges (pm_ij is the probability of an edge between two nodes belonging to i and j)
pm <- matrix (c(.3, .001, .001, .003,
                .001, .2, .005, .002,
                .001, .005, .2, .001,
                .003, .002, .001, .3), nrow=4, ncol=4)
g_sbm <- sample_sbm(100, pref.matrix=pm, block.sizes=c(25,25,25,25))
#g <- sample_sbm(200, pref.matrix=pm, block.sizes=c(50,50,50,50))
E(g_sbm)$weight <- 1

sample_wsbm_poisson <- function(pref.matrix, block.sizes){
    n <- sum(block.sizes)
    n_blocks <- length(block.sizes)
    block_index <- unlist( mapply(rep, 1:n_blocks, block.sizes, SIMPLIFY = FALSE) )
    adj_matrix <- matrix(0,ncol=n, nrow=n)
    for (i in 1:n){
        for (j in i:n){
            if (i==j) next
            adj_matrix[i,j]= rpois( 1, pref.matrix[block_index[i], block_index[j]] )
        }
    }
    graph_from_adjacency_matrix(adj_matrix, weighted=TRUE, mode="undirected")
}
g_wsbm_poisson <- sample_wsbm_poisson(5*pm,c(25,25,25,25))

sample_wsbm_binomial <- function(pref.matrix, block.sizes){
    #unfinished
    n <- sum(block.sizes)
    n_blocks <- dim(pref.matrix)[1]
    total_weight <- 0
    for (i in 1:n_blocks){
        for (j in 1:n_blocks){
            total_weight <- total_weight + block.sizes[i] * block.sizes[j] * pref.matrix[i,j] / 2
        }
    }
    
    block_index <- unlist( mapply(rep, 1:n_blocks, block.sizes, SIMPLIFY = FALSE) )
    prob_matrix <- matrix(0,ncol=n, nrow=n)
    for (i in 1:n){
        for (j in i:n){
            if (i==j) next
            prob_matrix[i,j]= pref.matrix[block_index[i], block_index[j]]
        }
    }
    prob_vector <- as.vector(prob_matrix)
    index_vector <- 1:(n*n) #indices to identify the edges when sampling
    s <- sample(n*n, size=total_weight, prob=prob_vector, replace=TRUE)
    weights_vector <- rep(0, n*n)
    for (i in s) weights_vector[i] <- weights_vector[i]+1
    adj_matrix <- matrix(weights_vector, nrow=n, ncol=n)
    graph_from_adjacency_matrix(adj_matrix, weighted=TRUE, mode="undirected")
}

multiply_diagonal <- function(lambda, M){
    for (i in 1:dim(M)[1])
        M[i,i] <- M[i,i] * lambda
    return(M)
}

parametrized_wsbm_binomial_list <- function(pref.matrix, block.sizes, params=10^seq(0,4)){
    pm_list <- lapply(params, multiply_diagonal, M=pref.matrix)
    g_list <- lapply(pm_list, sample_wsbm_binomial, block.sizes=block.sizes)
}

ground_truth_sbm <- function(blocks=c(40,25,25,10)){
    indices <- 1:length(blocks)
    unlist(mapply(rep, indices, blocks, SIMPLIFY=FALSE))
}

pm <- matrix (c(.03, .01, .01, .03,
                .01, .02, .05, .02,
                .01, .05, .02, .01,
                .03, .02, .01, .03), nrow=4, ncol=4)

g_wsbm_binomial <- sample_wsbm_binomial(multiply_diagonal(15, pm),c(40,25,25,10))
wsbm_gt_clustering <- ground_truth_sbm(c(40,25,25,10))
par(mar=c(0,0,0,0)+.05)
plot(g_wsbm_binomial, vertex.size=9, vertex.color=wsbm_gt_clustering, vertex.label=NA, 
     edge.width=(E(g_wsbm_binomial)$weight))


wsbm_parameters <- 10^seq(0,3,length.out=50)
block_sizes <- c(40, 25, 25, 10)
g_wsbm_p_list <- parametrized_wsbm_binomial_list(pm, c(40,25,25,10), 
                                                 params=wsbm_parameters)




