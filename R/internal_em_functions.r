# Author: Monica Golumbeanu <monica.golumbeanu@bsse.ethz.ch>

########################################
# Internal function init_kmeans: performs kmeans clustering on a given
# time-series matrix
# and estimates for each cluster a cubic smoothing spline using the ssanova
# function of the gss package.
# Used for initialization of the EM algorithm.
# input: data_table - matrix containing on each row the time series values
#        K - number of clusters
#        time_points - vector of time-points
# output: kmeans_clust - cluster allocation by kmeans
#         gss_obj_list - list object of class ssanova; contains the estimated
#         mixed effects model parameters;
########################
init_kmeans = function(data_table, K, time_points) {
    t = ncol(data_table)
    gss_obj_list = vector("list", K)
    # replace missing values by interpolation; only done for kmeans clustering.
    data_p = t(apply(data_table, MARGIN = 1, function(x) zoo::na.spline(c(x))))
    kmeans_clust = kmeans(data_p,K)$cluster

    # for each cluster, estimate the mean curve using gss ssanova function
    for (k in seq_len(K)) {
        cluster_genes = which(kmeans_clust == k)
        if (length(cluster_genes)>0) {
            cluster_data = data_table[cluster_genes,]
            gss_obj_list[[k]] = fit_ssanova(cluster_data, time_points)
        }
    }
    return(list(kmeans_clust = kmeans_clust, gss_obj_list = gss_obj_list))
}


##############################################
# Internal function calc_prob: given a data point and a model,
#                              calculates likelihood and posterior probabilities
# input: data_point: vector of size T where T is the number of time points
#        K: number of clusters
#        pi_k: mixing coefficients for each cluster component
#        gss_obj_list: list of estimated gss objects corresponding to a fitted
#                      smooting spline for each cluster; the length of the list
#                      is K.
# output: log_like: log-likelihood (proba of the data point given the model)
#         posterior: posterior probability of the model given the data point
#         clust: cluster assignment based on posterior probabilities
##############################################
calc_prob=function(data_point, K, pi_k, gss_obj_list) {
    t = length(data_point)
    components = vector("double", K)
    na_pos = is.na(data_point)
    for (k in seq_len(K)) {
        if(!is.null(gss_obj_list[[k]]$fit_model)) {
            # replacing the missing values with the mean of the normal
            # distribution; does not influence the result
            data_point[na_pos] = gss_obj_list[[k]]$est_mu[na_pos]
            # compute the covariance matrix Sigma
            Sigma_mat = matrix(data = 10^(gss_obj_list[[k]]$fit_model$zeta)*
                                    (gss_obj_list[[k]]$fit_model$varht),
                                    nrow = t, ncol = t) +
                        diag(x = gss_obj_list[[k]]$fit_model$varht,
                                nrow = t, ncol = t)

            components[k] = pi_k[k] *
                        dmvnorm(data_point, gss_obj_list[[k]]$est_mu, Sigma_mat)
        }
    }
    log_likelihood = log(sum(components))
    posterior = components/sum(components, na.rm = TRUE)

    # assign genes to clusters depending on the maximum posterior
    assigned_cluster = which(posterior == max(posterior))[1]

    return(list(log_like = log_likelihood,
                posterior = posterior, clust = assigned_cluster))
}

##############################################
# Internal function calc_ll: given a data point and a model,
#                              calculates the total log likelihood
# input: data_table: data frame with all the time-series data
#        K: number of clusters
#        pi_k: mixing coefficients for each cluster component
#        gss_obj_list: list of estimated gss objects corresponding to a fitted
#                      smooting spline for each cluster; K = length of the list
# output: log_like: log likelihood (proba of the data point given the model)
##############################################
calc_ll=function(data_table, K, pi_k, gss_obj_list) {
    vec_ll = vector("double", nrow(data_table))
    for (i in seq_len(nrow(data_table))) {
        components = vector("double", K)
        data_point = data_table[i,]
        na_pos = is.na(data_point)
        t = length(data_point)
        for (k in seq_len(K)) {
            if (length(gss_obj_list)>=k) {
                if(!is.null(gss_obj_list[[k]]$fit_model)) {
                    # replacing the missing values with the mean of the normal
                    # distribution; does not influence the result
                    data_point[na_pos] = gss_obj_list[[k]]$est_mu[na_pos]
                    # compute the covariance matrix Sigma
                    Sigma_mat = matrix(
                                data = 10^(gss_obj_list[[k]]$fit_model$zeta)*
                                        (gss_obj_list[[k]]$fit_model$varht),
                                            nrow = t, ncol = t) +
                                    diag(x = gss_obj_list[[k]]$fit_model$varht,
                                        nrow = t, ncol = t)
                    components[k] = pi_k[k] *
                        dmvnorm(data_point, gss_obj_list[[k]]$est_mu, Sigma_mat)
                }
            }
        }
        vec_ll[i] = log(sum(components))
    }
    ll = sum(vec_ll)
    return(ll)
}

##########################
# Internal function do_E_step: performs the expectation step of the EM
# input: data_table : matrix with the time-series
#        K: number of clusters
#        pi_k: mixing coefficients
#        gss_o_list: list of gss objects containing the estimated mixed-effects
#        models for each cluster
# output: mat_post: new posterior probabilities of data points to belong to
#           each cluster given model
#         cluster_assignment: vector with the assigned cluster for each point
##########################
do_E_step = function(data_table, K, pi_k, gss_obj_list) {
    cluster_assignment = vector("integer", K)
    N = nrow(data_table)
    mat_post = matrix(nrow = N, ncol = K)
    for (i in seq_len(N)) {
        data_point = data_table[i,]
        na_pos = is.na(data_point)
        t = length(data_point)
        components = vector("double", K)
        for (k in seq_len(K)) {
            if (length(gss_obj_list)>=k) {
                if(!is.null(gss_obj_list[[k]]$fit_model)) {
                    # replacing the missing values with the mean of the normal
                    # distribution; does not influence the result
                    data_point[na_pos] = gss_obj_list[[k]]$est_mu[na_pos]
                    # compute the covariance matrix Sigma
                    Sigma_mat = matrix(
                                data = 10^(gss_obj_list[[k]]$fit_model$zeta)*
                                        (gss_obj_list[[k]]$fit_model$varht),
                                            nrow = t, ncol = t) +
                                    diag(x = gss_obj_list[[k]]$fit_model$varht,
                                            nrow = t, ncol = t)
                    components[k] = pi_k[k] *
                        dmvnorm(data_point, gss_obj_list[[k]]$est_mu, Sigma_mat)
                }
            }
        }
        mat_post[i,] = components/sum(components, na.rm = TRUE)
        cluster_assignment[i] = which(mat_post[i,] == max(mat_post[i,]))[1]
    }

    return(list(mat_post = mat_post, cluster_assignment = cluster_assignment))
}


#########################################
# Internal function do_M_step: performs the maximisation step of the EM
# input: data_table: matrix with the time-series
#        time_points: vector of time-points
#        K: number of clusters
#        t_gss_obj_list: previously-estimated gss smoothing spline model
#        thresh: threshold to consider for the MC EM procedure
#        mat_post: posterior probability matrix
# output: new_gss_obj_list: new gss model solution of the M step
#########################################
do_M_step = function(data_table, time_points, K,
                        t_gss_obj_list, thresh, mat_post){
    # M-step
    new_gss_obj_list = vector("list", K)
    # number of elements and of time points
    N = nrow(data_table)
    t = ncol(data_table)
    for (k in seq_len(K)) {
        # use the rejection-controlled sampling to exclude low-probability
        # cluster members from
        # the gss smoothing spline estimation in order to obtain a more robust
        # estimation and avoid errors
        good_data_points = which(mat_post[,k]>=thresh)
        bad_data_points = which(mat_post[,k]<thresh)
        # with probability 1-p/thresh, assign posterior proba 0 to outliers
        reassignment = as.vector(sapply(mat_post[bad_data_points,k],
                                function(x, p = thresh){
                                    sample(c(1,0),1,prob=c(x/p,(1-x/p)))}))
        cluster_i = c(good_data_points, bad_data_points[reassignment == 1])
        if (length(cluster_i)>0){
            # obtain the weights for the penalized log-likelihood
            t_wei = t(as.matrix(mat_post[c(cluster_i),k]) %*%
                          matrix(1, nrow = 1, ncol = t))

            # maximize the penalized log-likelihood
            new_gss_obj_list[[k]] = fit_ssanova(data_table[cluster_i,],
                                                    time_points, t_wei)
        }
    }
    return(new_gss_obj_list)
}

##########################
# Internal function compute_pi_k: calculates mixing coefficients
# input: mat_post: posterior matrix
#        N: number of total observations
#        K: number of clusters
# output: pi_k
#########################
compute_pi_k = function(mat_post, N, K) {
    pi_k = vector("double", K)
    pi_k = colSums(mat_post, na.rm = TRUE)/N
    return(pi_k)
}

##########################
# Internal function do_EM : given a starting point, performs the EM iterations
# input: data_table: matrix containing the time-series data
#        time_points: vector with the time points
#        K: number of clusters
#        start_model: starting model configuration
#        em_iter_max: maximum number of iterations for EM
#        mc_em_iter_max: maximum number of iterations for the MC resampling
#        em_ll_convergence: convergence threshold for the likelihood difference
#                           between 2 consecutive iterations
# output: em_gss_obj_list: list of gss smoothing spline objects corresponding to
#                          estimated parameters for the smoothing splines
#         em_pi_k = estimated mixing coefficients
#         em_ll = vector containing the log likelihood after each EM iteration
##########################

do_EM = function(data_table, time_points, K, start_model, em_iter_max,
                    mc_em_iter_max, em_ll_convergence) {
    # number of elements and of time-points
    N = nrow(data_table)
    t = ncol(data_table)
    # initialize convergence check parameters and model parameters
    n_iter = previous_ll = 0
    converged = FALSE
    ll = em_pi_k = NULL
    thresh = 1
    decrease_iter = 0
    cluster_assignment = start_model$kmeans_clust
    em_gss_obj_list = start_model$gss_obj_list
    em_pi_k = as.vector(table(start_model$kmeans_clust))/N
    # perform initial E-step in order to get the initial posterior matrix
    E_step_list = do_E_step(data_table, K, em_pi_k, em_gss_obj_list)
    mat_post = E_step_list$mat_post

    #calculate likelihood of starting point
    n_iter = n_iter + 1
    previous_ll = calc_ll(data_table, K, em_pi_k, em_gss_obj_list)
    ll[n_iter] = previous_ll

    ## do EM iterations
    while(converged == FALSE & n_iter <= em_iter_max) {
        # E-step: compute new posterior matrix
        E_step_list = do_E_step(data_table, K, em_pi_k, em_gss_obj_list)
        # calculate mixing coefficients
        new_pi_k = compute_pi_k(E_step_list$mat_post, N, K)
        # MC-EM algorithm: subsample low-posterior genes while likelihood is
        # decreasing
        subsample = TRUE
        while(subsample) {
            # M-step: adjust the cluster weights and parameters
            new_gss_obj_list = do_M_step(data_table, time_points, K,
                                            em_gss_obj_list, thresh,
                                            E_step_list$mat_post)
            # evaluate the new log likelihood
            new_ll = calc_ll(data_table, K, new_pi_k, new_gss_obj_list)
            ll_dif = new_ll - previous_ll
            decrease_iter = decrease_iter + 1
            if (ll_dif >= - em_ll_convergence || decrease_iter > mc_em_iter_max)
            {
                subsample = FALSE
                message("log likelihood is increasing")
            }
        }
        if(ll_dif >= - em_ll_convergence) {
            n_iter = n_iter + 1
            ll[n_iter] = new_ll
            previous_ll = new_ll
            # update the parameters
            em_gss_obj_list = new_gss_obj_list
            mat_post = E_step_list$mat_post
            em_pi_k = new_pi_k
            cluster_assignment = E_step_list$cluster_assignment
            decrease_iter = 0
        }
        if ((abs(ll_dif) <= em_ll_convergence) ||
            (decrease_iter > mc_em_iter_max))  {
            converged = TRUE
        }
    }

    return(list(em_gss_obj_list = em_gss_obj_list, em_pi_k = em_pi_k,
                em_mat_post = mat_post,
                em_cluster_assignment = cluster_assignment,
                em_ll = ll, ts_data = data_table, ts_time_points = time_points))
}
