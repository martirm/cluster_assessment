#test cluster significance by applying clustering algorithms to the test networks and then evaluating
#scoring functions on them. For networks with a "ground-truth clustering", evaluate that one as well and compare

source("load_data.R")
library(clustAnalytics)
library(xtable)
library(pbapply)
library(mcclust)


set.seed(1)
significance_table_karate <- evaluate_significance_r(karate, alg_list=clust_alg_list,
                                                     weight_sel="max_weight", n_reps=100,
                                                     ground_truth=TRUE, gt_clustering=karate_gt_clustering,
                                                     w_max=NULL)
significance_table_forex <- evaluate_significance_r(g_forex, lower_bound=0, upper_bound=1, alg_list=clust_alg_list,
                                                    n_reps=100, w_max=1)
significance_table_news <- evaluate_significance_r(g_news, n_reps=10, w_max=NULL, alg_list=clust_alg_list)

significance_table_wsbm <- evaluate_significance_r(g_wsbm_binomial, ground_truth=TRUE,  alg_list=clust_alg_list,
                                                   gt_clustering=wsbm_gt_clustering, n_reps=100,
                                                   weight_sel="max_weight")
significance_table_enron <- evaluate_significance_r(weighted_enron, weight_sel="max_weight",  alg_list=clust_alg_list,
                                                    n_reps=100, w_max=NULL)

significance_table_sn <- evaluate_significance_r(g_social_network, weight_sel="max_weight",  alg_list=clust_alg_list,
                                                 n_reps=10, w_max=NULL)


####
significance_caption <- function(name) paste("Values of scoring functions for the", name,
                                             "graph (left), compared to those of their averages for 
                                             100 randomized samples (right), and its percentile rank
                                             (parentheses).")

print(xtable(add_arrow_rownames(significance_table_karate), caption = significance_caption("karate club"),
             label="significance_karate", align="lrrrrrr"), 
             latex.environments=c("center", "adjustbox"))
print(xtable(add_arrow_rownames(significance_table_forex), caption = significance_caption("Forex"), 
             label="significance_forex", align="lrrrrr"))
print(xtable(add_arrow_rownames(significance_table_wsbm), caption = significance_caption("weighted SBM"),
             label="significance_sbm", align="lrrrrrr"))
print(xtable(add_arrow_rownames(significance_table_news), caption = significance_caption("company news"),
             label="significance_news", align="lrrrrr"))
print(xtable(add_arrow_rownames(significance_table_sn), caption = significance_caption("social network"),
             label="significance_sn", align="lrrrrr"))
print(xtable(add_arrow_rownames(significance_table_enron), caption = significance_caption("enron"),
             label="significance_sn", align="lrrrrr"))
