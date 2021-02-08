source("stability_functions.R")
source("load_data.R")
library(clustAnalytics)
library(xtable)



b_karate <- boot_alg_list(alg_list, karate, return_data=TRUE) 
b_forex <- boot_alg_list(alg_list, g_forex, return_data=TRUE)
b_news <- boot_alg_list(alg_list, g_news, return_data=TRUE)
b_wsbm <- boot_alg_list(alg_list, g_wsbm_binomial, return_data=TRUE)
b_enron <- boot_alg_list(alg_list, weighted_enron, return_data=TRUE)


b_karate_cw <- boot_alg_list(alg_list, karate, return_data=TRUE, type="cluster-wise") 
b_forex_cw <- boot_alg_list(alg_list, g_forex, return_data=TRUE, type="cluster-wise")
b_news_cw <- boot_alg_list(alg_list, g_news, return_data=TRUE, type="cluster-wise")
b_wsbm_cw <- boot_alg_list(alg_list, g_wsbm_binomial, return_data=TRUE, type="cluster-wise")



b_rnd_karate <- boot_alg_list(alg_list, rewireCpp(karate, weight_sel="max_weight"), return_data=TRUE)
b_rnd_forex <- boot_alg_list(alg_list, rewireCpp(g_forex, upper_bound=1), return_data=TRUE)
b_rnd_news <- boot_alg_list(alg_list, rewireCpp(g_news, weight_sel="max_weight"), return_data=TRUE)
b_rnd_wsbm <- boot_alg_list(alg_list, rewireCpp(g_wsbm_binomial, weight_sel="max_weight"), return_data=TRUE)
b_rnd_enron <- boot_alg_list(alg_list, rewireCpp(weighted_enron, weight_sel="max_weight"), return_data=TRUE)


b_rnd_karate_cw <- boot_alg_list(alg_list, rewireCpp(karate, weight_sel="max_weight"), return_data=TRUE, type="cluster-wise")
b_rnd_forex_cw <- boot_alg_list(alg_list, rewireCpp(g_forex, upper_bound=1), return_data=TRUE, type="cluster-wise")
b_rnd_news_cw <- boot_alg_list(alg_list, rewireCpp(g_news_small, weight_sel="max_weight"), return_data=TRUE, type="cluster-wise")
b_rnd_wsbm_cw <- boot_alg_list(alg_list, rewireCpp(g_wsbm_binomial, weight_sel="max_weight"), return_data=TRUE, type="cluster-wise")







get_latex_table(b_karate, b_rnd_karate, "Zachary")
get_latex_table(b_forex, b_rnd_forex, "Forex")
get_latex_table(b_wsbm, b_rnd_wsbm, "wSBM")
get_latex_table(b_news, b_rnd_news, "News")
get_latex_table(b_enron, b_rnd_enron, "Enron")

jaccard_table(b_karate_cw, "Zachary")
jaccard_table(b_forex_cw, "Forex")
jaccard_table(b_wsbm_cw, "wSBM")
jaccard_table(b_news_cw, "News")
