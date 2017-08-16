library(TMixClust)
library(testthat)

# this test checks whether the clustering operation on the toy data works,
# in terms of how many genes are in each cluster.
test_that("test: clustering on the toy data gives the correct results", {
    t_obj = analyse_stability(toy_data_df, nb_clusters = 3,
                              nb_clustering_runs = 3, nb_cores = 1)
    c_sol = sum(sort(table(t_obj$em_cluster_assignment)) == c(30,30,31))
    expect_equal(c_sol, 3)
})


