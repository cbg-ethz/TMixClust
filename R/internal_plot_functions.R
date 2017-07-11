# Author: Monica Golumbeanu <monica.golumbeanu@bsse.ethz.ch>

###########################
# Internal function plot_log_like: plots the likelihood for all EM
# iterations to a file
# input: ll_vec: vector of log likelihood for all EM iterations
#        out_file: full path of the file where the plot will be saved
###########################
plot_log_like = function(ll_vec, out_file) {
    pdf(width=10, height=10, out_file, useDingbats=FALSE)
    par(mfcol=c(1,1))
    plot(ll_vec, type = "b", pch=20, lwd=2, xlab = "iteration",
            ylab = "log likelihood", cex.axis=1.3, cex.lab=1.5)
    dev.off()
}

###########################
# Internal function plot_clusters: generates nb_clusters plots, where
#                                  nb_clusters is the number of clusters;
#                                  each plot shows the time-series in the
#                                  respective cluster
# input: data_table: matrix containing the time-series
#        time_points: the values of the time-points
#        nb_clusters: number of clusters
#        final_clust: cluster assignment of the data points
#        out_file: full path of the .pdf file containing the plots
#        curve_file: file to be used for retrieving the average smoothing
#                       spline curves
#        data_color: color to be used for plotting the time-series
#        x_label: label to be used for the x axis
#        y_label: label to be used for the y axis
#        data_ticks: ticks to be used for the plots
#########################
plot_clusters = function(data_table, time_points, nb_clusters, final_clust,
                            report_folder, data_color, x_label, y_label) {
    data_color = adjustcolor(data_color, alpha.f = 0.5)
    yl = range(data_table, na.rm = TRUE) + c(-1,1)
    pdf(width=(5*nb_clusters), height=5,
            paste(report_folder,"clusters/cluster_plots.pdf",sep=""),
            useDingbats=FALSE)
    par(mfcol=c(1,nb_clusters),mar=c(5,5,5,5), oma=c(0,0,0,0))
    for (k in 1:nb_clusters) {
        genes_k = which(final_clust==k)
        for (j in 1:length(genes_k)) {
            plot(time_points, data_table[genes_k[j],], ylim=yl,
                    axes = FALSE, xlab="", ylab="", main="",
                    col=data_color, type="l", pch=20,lwd=2)
            par(new=TRUE)
        }
        # plot the fitted curve
        if (file.exists(paste(report_folder,
                            "estimated_curves/estimated_curve_cluster_",
                            k,".txt",sep=""))) {
            est_curve = read.table(paste(
                report_folder,"estimated_curves/estimated_curve_cluster_",
                k,".txt",sep=""))
            par(new=TRUE)
            plot(as.numeric(est_curve[4,]), as.numeric(est_curve[1,]),
                    ylim=yl, axes = FALSE, xlab="", ylab="", main="",
                    col="black", type="l", lwd=2, lty=2)
            par(new=TRUE)
            plot(as.numeric(est_curve[4,]), as.numeric(est_curve[2,]),
                    ylim=yl, axes = FALSE, xlab="", ylab="", main="",
                    col="black", type="l", lwd=2, lty=2)
            par(new=TRUE)
            plot(as.numeric(est_curve[4,]), as.numeric(est_curve[3,]), ylim=yl,
                    axes = FALSE, xlab="", ylab="", main="",
                    col="black", type="l", lwd=2)
            par(new=FALSE)
            # add axes to plot
            axis(side = 1, cex.axis=2.5,cex.lab=2.5) #at = data_ticks,
            axis(side = 2,cex.axis=2,cex.lab=2)
            title(main=paste("cluster ", k, sep=""), xlab = x_label,
                    ylab = y_label, cex.lab=2, cex.main=2)

        }
    }
    dev.off()
}



