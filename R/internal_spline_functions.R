# Author: Monica Golumbeanu <monica.golumbeanu@bsse.ethz.ch>

####################
# Internal function create_z - constructs the mixed effects unity matrix z
# needed by ssanova for fitting the spline
# Input: data_table - table containing the data, can have missing values
# Output: z - the constructed matrix
####################
create_z = function(data_table) {
    nb_numeric = apply(data_table, MARGIN = 1, function(x) sum(1-is.na(c(x))))
    nb_elements = sum(nb_numeric)
    previous = 0
    z = c(rep(1,nb_numeric[1]), rep(0,nb_elements-nb_numeric[1]))
    if(nrow(data_table)>1) {
        for (i in 2:nrow(data_table)) {
            previous = previous + nb_numeric[i-1]
            column = c(rep(0,previous), rep(1,nb_numeric[i]),
                        rep(0,nb_elements-previous-nb_numeric[i]))
            z = cbind(z,column)
        }
    }
    colnames(z) = NULL
    return(z)
}

#######################
# Function fit_ssanova
# Given a data table, uses the gss package to fit a cubic spline
# Input: data_to_fit - matrix of time-series to use for fitting the curve
#        time_points - list of time points
# Output: fit_model - the ssanova object containing the fitted model
#         est_mu - the estimated mean curve
#######################
# correlation matrix construction for mixed effects models, needed by ssanova
# the matrix is a square matrix of size equal the number of elements (env)
mixed_e_fun = function(zeta, env) {
    diag(10^(-zeta), env)
}

fit_ssanova = function(data_to_fit, time_points, wei) {
    # initialisations
    fit_model = est_mu = NULL
    data_size = nrow(data_to_fit)
    data_vector = as.numeric(as.vector(t(data_to_fit)))
    selected_numeric = which(as.numeric(is.na(data_vector))==0)
    time_grid = rep(time_points,data_size)
    tm = time_grid[selected_numeric]
    # specify the random effects part of the model
    random_sigma = list(fun = mixed_e_fun, env = data_size)
    random_effect = list(z = create_z(data_to_fit),
                            sigma = random_sigma, init = -1)
    # specify the fixed effects part of the model
    fixed_custom <- list(nphi=1, mkphi=mkphi.cubic, mkrk=mkrk.cubic,
                            env=c(min(time_points),max(time_points)))

    # fit the SSANOVA model with or without specified weights
    if(missing(wei)) {
        tryCatch({
            #print("fitting without weights")
            # fit model
            fit_model = ssanova(data_vector[selected_numeric]~tm,
                                    random=random_effect,
                                    type=list(tm = list("custom",fixed_custom)),
                                    alpha=1)
            # estimate the mean curves and model parameters
            est_mu = predict(fit_model,
                                data.frame(tm = time_points,
                                        offset=rep(0,length(time_points))))
        }, warning = function(w) {
            print(w)
        }, error = function(e) {
            print(e)
            stop("There was a problem with fitting the SSANOVA model. Try
running again or reduce the amount of missing values in the
data if applicable.")
        })
    } else {
        tryCatch({
            #print("fitting with weights")
            # specify the weights
            wei_vector = as.numeric(as.vector(wei))
            # fit the model
            fit_model = ssanova(data_vector[selected_numeric]~tm,
                                random=random_effect,
                                type=list(tm = list("custom",fixed_custom)),
                                weights=wei_vector[selected_numeric],alpha=1)
            # estimate the mean curves and model parameters
            est_mu = predict(fit_model, data.frame(tm = time_points,
                                            offset=rep(0,length(time_points))))
        }, warning = function(w) {
            print(w)
        }, error = function(e) {
            print(e)
            stop("There was an estimation problem with fitting the
SSANOVA model. Try running again or reduce the amount of missing values in
the data if applicable.")
        })
    }
    return(list(fit_model = fit_model, est_mu = est_mu))
}

