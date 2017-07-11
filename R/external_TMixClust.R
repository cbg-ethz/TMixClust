# Author: Monica Golumbeanu <monica.golumbeanu@bsse.ethz.ch>

##############################
#' @title Extracts a time series data frame from a text file
#'
#' @description \code{get_time_series_df} creates a data frame containing time
#' series data from a file.
#' @param data_file path to a tab-delimited text file containing the time series
#'  data formatted such that each row contains a time-series represented by its
#'  name (e.g. gene name, protein name, etc.) and the values at each time point.
#'
#' @return A data frame containing the time series
#'
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#'
#' @examples
#' # Load a simulated toy time-series data provided with the package
#' toy_data_file = system.file("extdata", "toy_time_series.txt",
#' package = "TMixClust")
#' toy_data= get_time_series_df(toy_data_file)
#'
#' # Print the first lines of the resulting data frame
#' print(head(toy_data))
#'
#' @export
#'
get_time_series_df = function(data_file) {
    if (!file.exists(data_file))
        stop("Input file does not exist.")
    if (missing(data_file))
        stop("A file containing the data to be loaded needs to be specified.")

    data_table = read.table(data_file, header=TRUE, na.strings =" ",
                            row.names=1, sep="\t")
    return(data_table)
}

##############################
#' @title Plots all the time series stored in a data frame object
#'
#' @description \code{plot_time_series_df} allows the user to visualise the time
#'  series from a given data set.
#'
#' @param ts_df data frame containing on each row a time-series
#' @param time_points vector containing the values of the time points.
#' Default: \code{c(1:ncol(time_series_df))}.
#' @param data_color color of the time series to be used for the plot.
#' Default is orange.
#' @param x_label label of the x axis of the plot. Default is "time"
#' @param y_label label of the y axis of the plot. Default is "value"
#' @param plot_title title of the plot. Default is "Time series plot".
#'
#' @return Plots a figure with all the the time series in the data set
#'
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#'
#' @examples
#' # Load the toy time series data provided with the TMixClust package
#' data(toy_data_df)
#'
#' # Plot the time series
#' plot_time_series_df(toy_data_df)
#' @export
#'
plot_time_series_df = function(ts_df, time_points = c(1:ncol(ts_df)),
                            data_color="#fd8d3c", x_label="time",
                            y_label="value", plot_title="Time series plot") {
    if (missing(ts_df))
        stop("A data-frame containing the time series needs to be provided.")

    # Plot the time series
    par(mar=c(5,5,5,5), oma=c(0,0,0,0))
    data_color = adjustcolor(data_color, alpha.f = 0.5)
    yl = range(ts_df, na.rm = TRUE) + c(-1,1)
    for (j in 1:nrow(ts_df)) {
        plot(time_points, ts_df[j,], ylim=yl, axes = FALSE, xlab="", ylab="",
            main="", col=data_color, type="l", pch=20,lwd=2)
        par(new=TRUE)
    }
    par(new=FALSE)
    # add axes to the plot
    axis(side = 1, cex.axis=2.5,cex.lab=2.5)
    axis(side = 2,cex.axis=2,cex.lab=2)
    title(main=plot_title, xlab = x_label, ylab = y_label,
            cex.lab=2, cex.main=2)
}

#############################
#' @title Clusters the time series data in a given number of groups
#'
#' @description \code{TMixClust} is the central function of the package.
#' It clusters the given time series data into a specified number of clusters.
#'
#' @param time_series_df data frame containing the time series.
#' Each row is a time series comprised of the time series name which is also the
#'  row name, and the time series values at each time point.
#' @param time_points vector containing numeric values for the time points.
#' Default: \code{c(1:ncol(time_series_df))}.
#' @param nb_clusters desired number of clusters
#' @param em_iter_max maximum number of iterations for the
#' expectation-maximization (EM) algorithm. Default: 1000.
#' @param mc_em_iter_max maximum number of iterations for Monte-Carlo
#' resampling. Default is 100.
#' @param em_ll_convergence convergence threshold for likelihood improvement.
#' Default is 0.001.
#'
#' @return list object with the following attributes:
#' \itemize{
#' \item \code{em_gss_obj_list} object of class \code{gss} containing
#' estimated parameters of the mixed-effects model
#' (see package vignette for more details).
#' \item \code{em_pi_k} vector containing the mixing coefficients
#' corresponding to each cluster
#' \item \code{em_mat_post} matrix containing the posterior values for each
#' time series and cluster
#' \item \code{em_cluster_assignment} vector with the clustering attribution
#' for each time series
#' \item \code{el_ll} vector containing the log likelihood values at each
#' iteration in the EM algorithm
#' \item \code{ts_data} the same as the input time series data-frame
#' \item \code{ts_time_points} the same as the input time-points vector
#' }
#'
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#'
#' @importFrom mvtnorm dmvnorm
#' @import gss
#' @import stats
#' @importFrom utils read.table
#'
#' @examples
#' # Load the toy time series data provided with the TMixClust package
#' data(toy_data_df)
#'
#' # Cluster the toy data with default parameters
#' TMixClust_obj = TMixClust(toy_data_df)
#'
#' @export
#'
TMixClust = function(time_series_df, time_points = c(1:ncol(time_series_df)),
                    nb_clusters = 2, em_iter_max = 1000, mc_em_iter_max = 10,
                    em_ll_convergence = 0.001) {
    # input checking conditions
    if (missing(time_series_df))
        stop("A data frame containing time series needs to be provided.")
    if(!length(time_points)==ncol(time_series_df))
        stop("The provided time series vector has different length than the
                number of columns in the time series data frame.")
    if(!(nb_clusters==round(nb_clusters) && nb_clusters > 1))
        stop("The number of clusters has to be a positive integer
                larger than 1.")
    if(!(em_iter_max==round(em_iter_max) && em_iter_max > 1))
        stop("The number of iterations for EM algorithm has to be a positive
                integer larger than 1.")
    if(!(mc_em_iter_max==round(mc_em_iter_max) && mc_em_iter_max > 1))
        stop("The number of iterations for the MC EM has to be a positive
                integer larger than 1.")
    if(em_ll_convergence < 0)
        stop("The convergence threshold for the likelihood convergence has to
                be positive.")

    # obtain starting configuration for EM with kmeans
    print("Initializing ...")
    k_init = init_kmeans(time_series_df, nb_clusters, time_points)

    # perform EM to estimate the parameters of the model and cluster the time
    # series based on their posterior
    print("Performing EM, this operation might take a few minutes ...")
    em_clust = do_EM(time_series_df, time_points, nb_clusters, k_init,
                    em_iter_max, mc_em_iter_max, em_ll_convergence)
    print("Clustering done, results stored in a TMixClust object.")

    return(em_clust)
}

#############################
#' @title Generates a series of files containing a summary of the TMixClust
#' analysis results
#'
#' @description \code{generate_TMixClust_report}
#'
#' @param TMixClust_object list object created by the \code{TMixClust}
#' function (see function \code{TMixClust})
#' @param report_folder full path of the folder where the report files will
#' be saved. Default is TMixClust_report/ folder in current working directory.
#' @param data_color color of the time series to be used when generating the
#' cluster plots. Default is orange.
#' @param x_label label of the x axis for the cluster plots. Default is "time"
#' @param y_label label of the y axis for the cluster plots. Default is "value"
#'
#' @return Produces a series of files containing information about the
#' clustering results and saves them in the provided folder location.
#' The folder contains the following:
#' \itemize{
#' \item \code{log-lihelihood.txt} - file with the log likelihood values at
#' each iteration on separate lines
#' \item \code{log-likelihood.pdf} - plot of log-likelihood at each iteration
#' \item \code{posterior.txt} - file with the posterior probabilities of all
#' the time-series for each cluster
#' \item \code{estimated_curves/} - folder containing a number of files equal to
#' the number of clusters; each file has 4 lines consisting of curve values and
#'  their confidence intervals (first 3 lines) for a discrete time grid
#'  (last line).
#' \item \code{clusters/} - folder containing a plot with the time series in
#' each cluster, a silhouette plot of the clustering configuration, as well as,
#' for each cluster, a file containing the names of the time series in the
#' respective cluster and a file containing the names and time series values
#' for the time series in each cluster.
#' }
#'
#' @examples
#' \dontrun{
#' # Load the toy time series data provided with the TMixClust package
#' data(toy_data_df)
#'
#' # Cluster the toy data with default parameters
#' TMixClust_obj = TMixClust(toy_data_df)
#'
#' # Generate a TMixClust report in the current working directory
#' generate_TMixClust_report(TMixClust_obj)
#' }
#'
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#'
#' @import gss
#' @import stats
#' @import cluster
#' @importFrom utils read.table write.table
#' @importFrom grDevices adjustcolor col2rgb dev.off pdf
#' @importFrom graphics abline axis par plot title
#'
#' @export
#'
generate_TMixClust_report = function(TMixClust_object, report_folder=paste(
    getwd(),"/TMixClust_report/",sep=""), data_color="#fd8d3c", x_label="time",
    y_label="value") {

    if (missing(TMixClust_object))
        stop("A TMixClust object needs to be provided.")
    if (is.null(TMixClust_object))
        stop("The TMixClust object is empty.")
    tryCatch( {
        col2rgb(data_color)
    } , error = function(e) {
        print(e)
        stop("Provided data color is not a valid color character.")
    })

    # extract number of clusters
    nb_clusters = length(TMixClust_object$em_pi_k)

    # check folder name
    if (!substring(report_folder, nchar(report_folder)) == "/") {
        report_folder = paste(report_folder,"/", sep="")
    }
    # create the necessary folders
    if (!dir.exists(report_folder)){
        tryCatch({
            print("creating report folder")
            dir.create(report_folder)
        }, warning = function(w) {
            print(w)
        }, error = function(e) {
            print(e)
            stop("There was a problem while creating the report folder.
                    Check if you have writing permission or the specified path
                    is correct.")
        })
    }
    if (!dir.exists(paste(report_folder,"clusters/", sep=""))){
        dir.create(paste(report_folder,"clusters/", sep=""))
    }
    if (!dir.exists(paste(report_folder,"estimated_curves/", sep=""))){
        dir.create(paste(report_folder,"estimated_curves/", sep=""))
    }

    # retrieve the curve with Bayesian confidence intervals
    for (k in 1:nb_clusters) {
        genes_k = which(TMixClust_object$em_cluster_assignment==k)
        if (length(genes_k)>0){
            if (!is.null(TMixClust_object$em_gss_obj_list[[k]]$fit_model)) {
                est_curve = matrix(0, nrow = 4, ncol = 100)
                time_grid = seq(min(TMixClust_object$ts_time_points),
                                max(TMixClust_object$ts_time_points),
                                length = 100)
                grid_predict = predict(
                    TMixClust_object$em_gss_obj_list[[k]]$fit_model,
                    data.frame(tm = time_grid), se.fit=TRUE)
                est_curve[1,] = grid_predict$fit + 1.96*grid_predict$se.fit
                est_curve[2,] = grid_predict$fit - 1.96*grid_predict$se.fit
                est_curve[3,] = grid_predict$fit
                est_curve[4,] = time_grid
                # write the curve and its confidence intervals to file
                write.table(est_curve,paste(
                    report_folder,"estimated_curves/estimated_curve_cluster_",k
                    ,".txt",sep=""), quote = FALSE, row.names = FALSE,
                    col.names = FALSE)
            }
        }
    }

    # write the cluster assignments to file and plot them
    for (k in 1:nb_clusters) {
        genes_k = which(TMixClust_object$em_cluster_assignment==k)
        write.table(TMixClust_object$ts_data[genes_k,],
                    paste(report_folder,"clusters/cluster",k,".txt",sep=""),
                    quote = FALSE)
        write.table(row.names(TMixClust_object$ts_data[genes_k,]),
                    paste(report_folder,"clusters/cluster",k
                            ,"_names.txt",sep=""), quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
    }
    plot_clusters(TMixClust_object$ts_data, TMixClust_object$ts_time_points,
                    nb_clusters,
                    TMixClust_object$em_cluster_assignment, report_folder,
                    data_color, x_label, y_label) #, data_ticks)

    # calculate silhouette coefficients for each cluster and
    # plot the distributions
    complete_ts = na.omit(TMixClust_object$ts_data)
    not_omitted = which(rownames(TMixClust_object$ts_data) %in%
                        rownames(complete_ts))
    m = data.matrix(complete_ts, rownames.force = FALSE)
    similarity_matrix = daisy(m)
    sil <- silhouette(TMixClust_object$em_cluster_assignment[not_omitted],
                        similarity_matrix)

    pdf(paste(report_folder,"clusters/silhouette.pdf",sep=""),
        useDingbats = FALSE)
    plot(sil, col = "#bdbdbd", border = "#bdbdbd",
            main = paste("Silhouette plot for K=",
                    nb_clusters, " clusters", sep=""))
    abline(v = mean(sil[,"sil_width"]), col="black", lwd=3, lty=2)
    dev.off()

    # write the likelihood values at each EM iteration to a file and plot it
    write.table(TMixClust_object$em_ll,
                paste(report_folder,"log-likelihood.txt",sep=""),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    plot_log_like(TMixClust_object$em_ll,
                    paste(report_folder,"log-likelihood.pdf",sep=""))

    # write the posterior probabilities to a file
    final_posterior = TMixClust_object$em_mat_post
    rownames(final_posterior) = rownames(TMixClust_object$ts_data)
    write.table(final_posterior, paste(report_folder,"posterior.txt",sep=""),
                    sep= "\t", quote = FALSE, col.names = FALSE)
}

#############################
#' @title Generates a silhouette plot for e given clustering configuration.
#'
#' @description \code{plot_silhouette}
#'
#' @param TMixClust_object list object created by the \code{TMixClust}
#' function (see function \code{TMixClust})
#' @param sim_metric character string taking one of the possible values:
#' "euclidean", "gower" or "manhattan". Default is "euclidean".
#' @param sil_color color of the bars representing the silhouette widths
#' on the plot
#'
#' @return List object with the following components:
#' \itemize{
#' \item \code{similarity_m} similarity matrix
#' \item \code{silh} silhouette object
#' }
#'
#' Renders a plot comprised of a set of barplots with the distributions of
#' silhouette
#' coefficients for the data points in each cluster. Each barplot has indicated
#' on its right hand side the total number of points in the corresponding
#' cluster.
#' The plot also indicates with a dotted line,
#' the overall average silhouette width,
#' whose value is specified at the bottom of the plot.
#'
#' @examples
#' # Load the TMixClust object associated to the toy time series data
#' # provided with the TMixClust package
#' data(best_clust_toy_obj)
#'
#' # Plot the silhouette for the clustering stored in the toy TMixClust object
#' plot_silhouette(best_clust_toy_obj)
#'
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#'
#' @import gss
#' @import stats
#' @import cluster
#' @importFrom grDevices adjustcolor col2rgb dev.off pdf
#' @importFrom graphics abline axis par plot title
#'
#' @export
#'
plot_silhouette = function(TMixClust_object, sim_metric="euclidean",
                            sil_color="#bdbdbd") {
    # extract the number of clusters
    nb_clusters = length(TMixClust_object$em_pi_k)

    # calculate silhouette coefficients for each cluster and
    # plot the distributions
    # this works only for the time-series without missing values,
    # therefore we have to first remove the data points with missing values:
    complete_ts = na.omit(TMixClust_object$ts_data)
    not_omitted = which(rownames(TMixClust_object$ts_data) %in%
                        rownames(complete_ts))
    m = data.matrix(complete_ts, rownames.force = FALSE)
    # compute similarity matrix and silhouette coefficients
    similarity_matrix = daisy(m, metric = sim_metric)
    sil <- silhouette(TMixClust_object$em_cluster_assignment[not_omitted],
                        similarity_matrix)
    # create silhouette plot
    plot(sil, col = sil_color, border = sil_color,
            main = paste("Silhouette plot for K=", nb_clusters,
                            " clusters", sep=""))
    abline(v = mean(sil[,"sil_width"]), col="black", lwd=3, lty=2)

    return(list(similarity_m = m, silh = sil))
}

#############################
#' @title  Stability analysis, clustering evaluation and optimal solution
#' selection
#'
#' @description \code{analyse_stability} Performs multiple clustering runs
#' with TMixClust, analyses the agreement between runs
#' with the Rand index and returns the clustering solution with the largest
#' likelihood.
#' A plot of agreement probability between all the runs and the run with the
#' maximum likelihood is produced.
#'
#' @param time_series_df data frame containing the time series.
#' Each row is a time series comprised of the time series name which is also
#' the row name, and the time series values at each time point.
#' @param time_points vector containing numeric values for the time points.
#' Default: \code{c(1:ncol(time_series_df))}.
#' @param nb_clusters desired number of clusters
#' @param em_iter_max maximum number of iterations for the
#' expectation-maximization (EM) algorithm. Default: 1000.
#' @param mc_em_iter_max maximum number of iterations for Monte-Carlo
#' resampling. Default is 100.
#' @param em_ll_convergence convergence threshold for likelihood improvement.
#' Default is 0.001.
#' @param nb_clustering_runs number of times the clustering procedure is
#' repeated on the input data. Default is 3.
#' @param nb_cores number of cores to be used to run the separate clustering
#' operations in parallel. Default is 2.
#'
#' @return TMixClust object with the highest likelihood.
#' Renders a plot showing the overall distribution of the Rand index, which
#' allows the user to assess clustering stability.
#'
#' @examples
#' # Load the toy time series data provided with the TMixClust package
#' data(toy_data_df)
#'
#' # Identify the most optimal clustering solution with 3 clusters
#' best_clust_obj = analyse_stability(toy_data_df, nb_clusters = 3,
#'                                    nb_clustering_runs = 4, nb_cores = 1)
#'
#' # Plot the time series from each cluster
#' for (i in 1:3) {
#'     # Extract the time series in the current cluster and plot them
#'     c_df=toy_data_df[which(best_clust_obj$em_cluster_assignment==i),]
#'     plot_time_series_df(c_df, plot_title = paste("cluster",i))
#' }
#'
#' @author Monica Golumbeanu, \email{monica.golumbeanu@bsse.ethz.ch}
#' @references Golumbeanu M, Desfarges S, Hernandez C, Quadroni M, Rato S,
#' Mohammadi P, Telenti A, Beerenwinkel N, Ciuffi A. (2017) Dynamics of
#' Proteo-Transcriptomic Response to HIV-1 Infection.
#'
#' @import doParallel
#' @import foreach
#' @importFrom flexclust randIndex
#' @importFrom grDevices adjustcolor col2rgb dev.off pdf
#' @importFrom graphics abline axis par plot title hist
#'
#' @export
#'
analyse_stability = function(time_series_df,
                                time_points = c(1:ncol(time_series_df)),
                                nb_clusters = 2, em_iter_max = 1000,
                                mc_em_iter_max = 100, em_ll_convergence = 0.001,
                                nb_clustering_runs=3, nb_cores=2){
    # input checking conditions
    if (missing(time_series_df))
        stop("A data frame containing time series needs to be provided.")
    if(!length(time_points)==ncol(time_series_df))
        stop("The provided time series vector has different length than
                the number of columns in the time series data frame.")
    if(!(nb_clusters==round(nb_clusters) && nb_clusters > 1))
        stop("The number of clusters has to be a positive integer larger
                than 1.")
    if(!(em_iter_max==round(em_iter_max) && em_iter_max > 1))
        stop("The number of iterations for EM algorithm has to be a
                positive integer larger than 1.")
    if(!(mc_em_iter_max==round(mc_em_iter_max) && mc_em_iter_max > 1))
        stop("The number of iterations for the MC EM has to be a positive
                integer larger than 1.")
    if(em_ll_convergence < 0)
        stop("The convergence threshold for the likelihood convergence has to
                be positive.")
    if(!(nb_clustering_runs==round(nb_clustering_runs) &&
        nb_clustering_runs >= 1))
        stop("The number of clustering runs has to be a positive integer
                strictly larger than 0.")
    if(!(nb_cores==round(nb_cores) && nb_cores >= 1))
        stop("The number of computing cores has to be a positive integer
                strictly larger than 0.")

    # set up the parallel environment
    registerDoParallel(cores=nb_cores)
    # perform each TMixClust clustering in parallel and store all the
    # TMixClust objects in a list
    cl_obj = foreach(i=1:nb_clustering_runs) %dopar%
            TMixClust(time_series_df, time_points, nb_clusters,
                        em_iter_max, mc_em_iter_max, em_ll_convergence)

    # find the solution with the highest likelihood
    max_pos = which.max(lapply(cl_obj, function(x) x$em_ll[length(x$em_ll)] ))
    best_clust_obj = cl_obj[[max_pos]]

    # calculate the rand index between each TMixClust object and the object with
    # the highest likelihood
    rand_ind = NULL
    for (i in 1:(nb_clustering_runs)) {
        g1 = cl_obj[[i]]$em_cluster_assignment
        g2 = cl_obj[[max_pos]]$em_cluster_assignment
        tab <- table(g1, g2)
        rand_ind[i] = randIndex(tab)
    }

    # plot the distribution of the obtained Rand indexes
    hist(rand_ind, xlab = "rand index",
            main = "Distribution of agreement probability",
            col = "#2c7fb8", breaks = 10)

    return(best_clust_obj)
}

