
# Author: Luis Cerdán <lcerdanphd_at_gmail.com>
#                     <https://github.com/lucerlab>

#######################################################################
#######################################################################
#####                                                             #####
#####       AUXILIARY FUNCTIONS FOR AUTO-D AND NEA-GMM            #####
#####                                                             #####
#######################################################################
#######################################################################


##########################################################################################################
# Function to determine outliers with Mahalanobis distance and the False Discovery Rate (FDR)
##########################################################################################################

find_outliers <- function(df,prob){
  ## Input:
  #   df: data frame containing VEE and f for all transitions (columns) and geometries (rows)
  #   prob: proportion of false outlier detections
  ## Output:
  #  df_clean: data frame without outliers
  #  outlier_geoms: list of outlier geometries  
  
  # copy original dataframe for later use
  df_tr <- df
  
  # store number of rows (n) and columns (p)
  n <- nrow(df)
  p <- ncol(df)

  # "normalize" distribution (take square root of f divided by VEE)
  df[,(p/2+1):p] <- sqrt(df[,(p/2+1):p]/df[,1:(p/2)])

  # standardize data (robustly)
  df <- apply(df,2,function(x){x-median(x)})
  df <- apply(df,2,function(x){x/IQR(x)})
  
  # compute sample median (mu.rob) and robust covariance matrix (S.rob)
  mu.rob <- robustbase::colMedians(df)

  M <-apply(df,2,function(x){x-median(x)}) # difference matrix
  S.rob <- (n-1)^{-1} * t(M) %*% M # covariance matrix
  
  # compute squared Mahalanobis distance
  D2 <- stats::mahalanobis(df,mu.rob,S.rob)
  
  # Obtain the p-values assuming Gaussianity (D2~chi^2_p)
  p.values <- 1 - pchisq(D2,p)
  
  # Sort them in increasing order
  sort.p.values <- sort(p.values,index.return=TRUE)
  
  # Which are the outliers? p.value < rho*prob/n
  sort.i_outliers <- which(sort.p.values$x < ((1:n)/n*prob))
  i_outliers <- sort.p.values$ix[sort.i_outliers]
  
  # store geometries containing outliers and remove them
  if(length(i_outliers)!=0){
    outlier_geoms <- rownames(df_tr)[i_outliers]
    df_clean <- df_tr[-c(i_outliers),]
  }else{
    df_clean <- df_tr
    outlier_geoms <- 'None'
  }
  
  # return clean dataframe and list of outliers
  return(list(df_clean,outlier_geoms))
}


##########################################################################################################
#######  Functions to compute absorption cross section for each transition and full spectrum
#######  from VEE and f pairs using auto-d
##########################################################################################################

# function to compute weighted Interquartile Range (wIQR)
weighted.IQR <- function (x, w) {
  ## Input:
  # x: data (Vertical transition energies)
  # w: weights
  ## Output:
  # wIQR: weighted IQR
  
  # normalize weights
  w <- w/sum(w)
  
  # Ensure x and w are in ascending order of x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  # find indices of quantiles 25% and 75%
  q1 <- which(cumsum(w)>0.25)[1]-1
  q2 <- which(cumsum(w)>0.75)[1]-1
  
  # find probabilities
  p1 <- 0.25-sum(w[1:q1])
  p2 <- 0.75-sum(w[1:q2])
  
  # Interpolate quantiles
  Q1 <- x[q1] + p1*(x[q1+1]-x[q1])
  Q3 <- x[q2] + p2*(x[q2+1]-x[q2])
  
  # Compute and return weighted IQR
  IQR_w <- Q3-Q1
  return(IQR_w)
}

# function to compute kernel width (eq. 8 of publication)
auto_d <- function(VEE,f){
  ## Input:
  # VEE: array with vertical energies
  # f: array with oscillator strengths
  ## Output:
  # d: optimized empirical bandwidth
  
  # number of geometries (structural configurations)
  confs <- length(VEE)
  
  # compute weighted standard deviation
  w_VEE <- VEE*f/confs # weight
  w_VEE <- w_VEE/sum(w_VEE) # normalized weight
  w_mu <- sum(w_VEE*VEE) # weighted mean
  w_var <- sum(w_VEE*(VEE-w_mu)^2) # weighted standard deviation
  w_IQR <- weighted.IQR(VEE,w_VEE) # weighted Interquartile Range
  A <- min(sqrt(w_var),w_IQR/1.34) # minimum of w_var and re-scaled IQR
  
  # compute and return d
  d <- 1.06*A*confs^(-1/5)
  return(d)
}


# function to compute spectrum (eq. 1 in publication) and confidence interval for one transition
band_gen <- function(VEE, f, E, ci){
  ## Input
  # VEE: array with vertical energies
  # f: array with oscillator strengths
  # E: array with energies to compute spectrum
  # ci: If true, computes confidence intervals
  ## Output
  # abs_cross: band cross section spectrum
  # abs_cross_ci: confidence interval
  # d: optimized empirical bandwidth
  
  # set physical constants and scaling factor for spectrum
  hbarev=6.582119514e-16  # Planck constant in eV * s
  c=299792458.0e0	        # Speed of light in vacuum in m/s
  eps0=8.854187817e-12    # Vacuum permitivity in C V-1 m-1
  me=9.10938356e-31       # electron mass in kg
  ec=1.6021766208e-19     # electron charge in C
  
  # store number of geometries (structural configurations)
  confs <- length(VEE)
  
  # count geometries with f = 0
  zeros <- which(f==0)
  
  # automatic selection of bandwidth
  if(length(zeros)>(confs-1)){ # if all geometries for a transition have f=0, fix d
    d <- 0.2}
  else{
    d <- 2*auto_d(VEE,f) # multiply by 2 because auto_d returns the half-width
  }
  
  # scaling factor
  factor_abs = 10000*pi*ec^2*hbarev/(2.0*me*c*eps0) # guaranties units of cm2 for cross section
  factor <- factor_abs/confs/d*sqrt(2/pi) # Gaussian scaling factor
  
  # compute spectrum
  abs_cross <- rep(0,length(E)) # array to store spectrum
  for(ii in 1:confs){
    abs_cross <- abs_cross + factor*f[ii]*VEE[ii]/E*exp(-2*(E-VEE[ii])^2/d^2)
  }
  
  # bootstrap to find band CI
  if(ci){
    n_boot <- 999
    spec_boot <- matrix(0, nrow = n_boot, ncol = length(E)) # matrix to store bootstrapped bands
    for(bb in 1:n_boot){
      id_boot <- sample(confs, replace = TRUE) # randomly sample conf observations with replacement
      f_boot <- f[id_boot]
      VEE_boot <- VEE[id_boot]
      for(ii in 1:confs){ # compute spectrum
        spec_boot[bb,] <- spec_boot[bb,] + factor*f_boot[ii]*VEE_boot[ii]/E*exp(-2*(E-VEE_boot[ii])^2/d^2)
      }
    }
    
    # compute lower and upper bound errors (95% CI) by taking quantiles (0.025, 0.975)
    # and subtracting original mean
    abs_cross_lci <- matrixStats::colQuantiles(spec_boot, probs = 0.025)
    abs_cross_lci <- abs(abs_cross_lci - abs_cross)
    abs_cross_uci <- matrixStats::colQuantiles(spec_boot, probs = 0.975)
    abs_cross_uci <- abs(abs_cross_uci - abs_cross)
  } else{
    abs_cross_lci <- rep(0,length(E))
    abs_cross_uci <- rep(0,length(E))
  }
  
  # return band spectrum, band errors, and empirical d
  return(list(abs_cross,abs_cross_lci,abs_cross_uci,d))
}


# function to compute auto-d spectrum
auto_d_spectra <- function(df,E,ci){
  ## Input:
  # df: data frame containing VEE and f for all transitions (columns) and geometries (rows)
  # E: array with energies to compute spectrum
  # ci: If true, computes confidence intervals
  ## Ouput:
  # bands: dataframe with band cross section spectrum
  # bands_lci: dataframe with band cross section spectrum lower confidence interval
  # bands_uci: dataframe with band cross section spectrum upper confidence interval
  # abs_cross: array with full cross section spectrum
  # abs_cross_lci: array with full cross section spectrum lower confidence interval
  # abs_cross_uci: array with full cross section spectrum upper confidence interval
  # d: array with optimized empirical linewidths
  
  # store number of transitions/bands/states and confs
  n_states <- ncol(df)/2
  confs <- nrow(df)
  
  # compute all bands cross section spectrum (parallelized)
  band_par <- foreach(
    ii = 1:n_states,
    .combine = 'cbind'
  ) %dopar% {
    band_tmp <- band_gen(df[,ii],df[,(ii+n_states)], E, ci)
    # merge to avoid list of lists
    c(band_tmp[[1]],band_tmp[[2]],band_tmp[[3]],band_tmp[[4]])
  }
  
  # convert to matrix form (to allow for one transition computations)
  band_par <- as.matrix(band_par)
  
  # separate bands from bandwidth
  n_rows <- nrow(band_par)
  
  bands <- band_par[1:length(E),] # band cross section spectrum
  bands_lci <- band_par[(length(E)+1):(2*length(E)),] # band cross section spectrum CI
  bands_uci <- band_par[(2*length(E)+1):(3*length(E)),] # band cross section spectrum CI
  d_opt <- band_par[n_rows,] # array with auto-d linewidths
  
  # generate full spectrum by adding up all bands contributions
  if(ncol(band_par)!=1){ # if there are multiple bands
    abs_cross <- rowSums(bands) # eq. 3 in publication
    abs_cross_lci <- sqrt(rowSums(bands_lci^2)) # eq. 5 in publication
    abs_cross_uci <- sqrt(rowSums(bands_uci^2)) # eq. 5 in publication
  } else {# if there is only one
    abs_cross <- bands
    abs_cross_lci <- bands_lci
    abs_cross_uci <- bands_uci
  }
  
  # return bands spectra and errors, full spectrum and errors, and linewidths
  return(list(bands, bands_lci, bands_uci, abs_cross, abs_cross_lci, abs_cross_uci, d_opt))
}


##########################################################################################################
#######  Functions to compute absorption cross section for each transition and full spectrum
#######  from VEE and f pairs using GMM-NEA
##########################################################################################################

# function to compute values of Normal distribution PDF at a given energy
gauss <- function(E,mu,var){
  ## Input:
  # E: energy at which to evaluate PDF
  # mu: array with means
  # var: array with variances
  ## Output:
  # g: array of PDF values for each mu and var
  E <- rep(E,length(mu))
  
  # compute and return g
  g <- 1/(sqrt(2*pi*var))*exp(-(E-mu)^2/(2*var))
  return(g)
}

# function to compute spectrum from GMM parameters (eqs 21-24 of publication)
gmm_band <- function(gmm_model, Theta_0, E, mu_vee, sd_vee, mu_m, sd_m){
  ## Input:
  # gmm_model: fitted mclust model
  # Theta_0: proportion of geometries with M=0
  # E: array with energies to compute spectrum
  # mu_vee: vector of means of input VEE values
  # sd_vee: vector of standard deviations of input VEE values
  # mu_m: vector of means of input M values
  # sd_m: vector of standard deviations of input M values
  ## Ouput:
  # abs_cross: band spectrum
  
  # set physical constants and scaling factor for spectrum
  hbarev=6.582119514e-16
  ec=1.6021766208e-19
  me=9.10938356e-31
  c=299792458.0
  eps0=8.854187817e-12
  J2eV=(3e8)^4*1.78266192e-36*1e42 # J to eV conversion factor
  factor_abs <- 10000*pi/(3*hbarev*c*eps0)/J2eV # guaranties units of cm2 for cross section
  
  # store gmm parameters and rescale to have original scales
  pi_k <- gmm_model$parameters$pro # weights/priors/coefficients 
  mu_e <- gmm_model$parameters$mean[1,]*sd_vee+mu_vee # means (VEE)
  mu_m <- gmm_model$parameters$mean[2,]*sd_m+mu_m # means (M)
  var_e <- gmm_model$parameters$variance$sigma[1,1,]*sd_vee^2 # variances (VEE)
  var_m <- gmm_model$parameters$variance$sigma[2,2,]*sd_m^2 # variances (M)
  sigma_em <- gmm_model$parameters$variance$sigma[1,2,]*sd_vee*sd_m # covariances, sigma_em = rho sigma1 sigma2
  rho <- sigma_em/sqrt(var_e*var_m) # correlation coefficients, rho = sigma_em /sqrt(sigma_e^2*sigma_m^2)
  # rescale mixture weights to account for zeros
  pi_k <- pi_k*(1-Theta_0)
  
  # compute band cross section
  abs_cross <- rep(0,length(E))
  for(ii in 1:length(E)){
    # replicate energy value for each mixture
    x1 <- rep(E[ii], length(pi_k))
    
    # compute mean and variance of conditional PDF
    mu_xy <- mu_m+rho*sqrt(var_m/var_e)*(x1-mu_e)
    var_xy <- (1-rho^2)*var_m
    
    # compute second order moment of M
    mean_fxy <- mu_xy^2 + var_xy
    
    # get PDF value at E
    phi <- gauss(E[ii],mu_e,var_e)
    
    # compute band shape
    abs_cross[ii] <- E[ii]*sum(pi_k*mean_fxy*phi)
  }
  
  # rescale spectrum
  abs_cross <- factor_abs*abs_cross
  
  return(abs_cross)
}

# function to fit GMM model and compute spectrum and confidence interval for one transition
band_gen_pdf <- function(df, E, ci){
  ## Input:
  # df: data frame containing VEE and M for given transition (columns) and geometries (rows)
  # E: array with energies to compute spectrum
  # ci: If true, computes confidence intervals
  ## Ouput:
  # abs_cross: band spectrum
  # abs_cross_lci: band spectrum lower confidence interval
  # abs_cross_uci: band spectrum upper confidence interval
  # k_opt: optimal number of mixtures
  # model_opt: model name (constraints)
  
  # store number of geometries (structural configurations)
  n_geoms <- dim(df)[1]
  
  # remove geometries with M = 0
  zeros <- which(df[,2]==0)
  #zeros <- c()
  
  if(length(zeros)>(n_geoms-5)){ # if there is less than 5 values with M != 0, set abs to 0
    abs_cross <- rep(0,length(E))
    abs_cross_lci <- rep(0,length(E))
    abs_cross_uci <- rep(0,length(E))
    k_opt <- 0
    model_opt <- 'none'
  }
  else{
    if(length(zeros)>0){ # if there are geometries with M=0, remove them
      df <- df[-zeros, ]
    }
    
    # store probability of being zero
    Theta_0 <- length(zeros)/n_geoms
    
    # Determine maximum number of clusters (k_max) to use:
    # Number of parameters of 2D GMM with full covariance matrix: n=K*Df-1, with Df <- (D*D - D)/2 + 2*D + 1; D=2
    # maximum K: k_max=(n-1)/Df/5. Use at least 5 observations per parameter to avoid overfitting.
    D <- 2 # dimensionality of data 
    Df <- (D*D - D)/2 + 2*D + 1 # Degrees of freedom for single normal
    k_max <- floor((dim(df)[1]+1)/Df/5) # use at least 5 observations per parameter
    k_max <- max(c(2,k_max)) # set minimum to 2 to force more flexibility
    k_max <- min(c(10,k_max)) # set maximum to 10 to avoid lengthy computations
    k_min <- 2 # set minimum to force more flexibility
    
    # standardize data (it is convenient to avoid excessive differences in VEE and M scales)
    df_stand <- df
    mu_vee <- mean(df[,1]) # mean
    sd_vee <- sd(df[,1]) # standard deviation
    df_stand[,1] <- (df[,1]-mu_vee)/sd_vee
    mu_m <- mean(df[,2]) # mean
    sd_m <- sd(df[,2]) # standard deviation
    df_stand[,2] <- (df[,2]-mu_m)/sd_m
    
    # Run Mclust model selection and fit&store optimal model and number of mixtures
    nclust <- mclustBIC(df_stand, G=k_min:k_max, verbose = FALSE)
    single_gmm <- Mclust(df_stand, x = nclust, verbose = FALSE)
    model_opt <- single_gmm$modelName
    k_opt <-single_gmm$G
    
    # compute transition cross section
    abs_cross <- gmm_band(single_gmm, Theta_0, E, mu_vee, sd_vee, mu_m, sd_m)
    
    # bootstrap to find CI
    if(ci){
      n_boot <- 999
      abs_boot <- matrix(NA, nrow = n_boot, ncol = length(E)) # matrix to store bootstrapped bands
      for(bb in 1:n_boot){
        id_boot <- sample(nrow(df), replace = TRUE) # randomly sample conf observations with replacement
        df_boot <- df[id_boot,]
        
        # restandardize data
        mu_vee <- mean(df_boot[,1]) # mean
        sd_vee <- sd(df_boot[,1]) # standard deviation
        df_boot[,1] <- (df_boot[,1]-mu_vee)/sd_vee
        mu_m <- mean(df_boot[,2]) # mean
        sd_m <- sd(df_boot[,2]) # standard deviation
        df_boot[,2] <- (df_boot[,2]-mu_m)/sd_m
        
        # refit mclust model
        # add tryCatch to catch exceptions (for some bootstrap replicas, mclust may fail)
        skip_to_next <- FALSE
        tryCatch(single_gmm <- Mclust(df_boot, x = nclust, verbose = FALSE), error = function(e){skip_to_next <- TRUE})
        if(skip_to_next){ next }
        
        # compute transition cross section
        abs_boot[bb,] <- gmm_band(single_gmm, Theta_0, E, mu_vee, sd_vee, mu_m, sd_m)
      }
      
      # compute lower and upper bound errors (95% CI) by taking quantiles (0.025, 0.975)
      # and subtracting original mean
      abs_cross_lci <- matrixStats::colQuantiles(abs_boot, probs = 0.025)
      abs_cross_lci <- abs(abs_cross_lci - abs_cross)
      abs_cross_uci <- matrixStats::colQuantiles(abs_boot, probs = 0.975)
      abs_cross_uci <- abs(abs_cross_uci - abs_cross)
    } else{
      abs_cross_lci <- rep(0,length(E))
      abs_cross_uci <- rep(0,length(E))
    }
  }
  
  # return band shape, band shape CI, optimized number of mixtures, and optimized model name
  return(list(abs_cross,abs_cross_lci,abs_cross_uci,k_opt,model_opt))
}

# function to compute GMM-NEA spectrum
gmm_nea_spectra <- function(df,E,ci){
  ## Input:
  # df: data frame containing VEE and f for all transitions (columns) and geometries (rows)
  # E: array with energies to compute spectrum
  # ci: If true, computes confidence intervals
  ## Ouput:
  # bands: dataframe with band cross section spectrum
  # bands_lci: dataframe with band cross section spectrum lower confidence interval
  # bands_uci: dataframe with band cross section spectrum upper confidence interval
  # abs_cross: array with full cross section spectrum
  # abs_cross_lci: array with full cross section spectrum lower confidence interval
  # abs_cross_uci: array with full cross section spectrum upper confidence interval
  # k_opt: optimal number of mixtures
  # model_opt: model name (constraints)
  
  # store number of df columns
  n_col <- ncol(df)
  
  # store number of transitions/bands/states
  n_states <- ncol(df)/2
  
  # Transform f to M (in Debye)
  # M = sqrt(J2eV*3*hbarev^2*ec^2/(2*me) * f/E)
  # J2eV=(3e8)^4*1.78266192e-36*1e42
  hbarev=6.582119514e-16 # eV s
  ec=1.6021766208e-19 # C
  me=9.10938356e-31 # kg
  J2eV=(3e8)^4*1.78266192e-36*1e42 # J to eV conversion factor
  factor_s = J2eV*3*hbarev**2*ec^2/(2*me)
  df[,(n_states+1):n_col] <- sqrt(factor_s*df[,(n_states+1):n_col]/df[,1:(n_states)])
  
  # compute all bands cross section spectrum (parallelized)
  band_par <- foreach(
    ii = 1:n_states, 
    .combine = 'cbind',
    .packages='mclust'
  ) %dopar% {
    band_tmp <- band_gen_pdf(df[,c(ii,ii+n_states)], E, ci)
    # merge to avoid list of lists
    c(band_tmp[[1]],band_tmp[[2]],band_tmp[[3]],band_tmp[[4]],band_tmp[[5]])
  }
  
  # convert to matrix form (to allow for one transition computations)
  band_par <- as.matrix(band_par)
  
  # separate bands from model name and convert to numeric values
  n_rows <- nrow(band_par)
  
  bands_plus_k <- band_par[1:(n_rows-1),]
  class(bands_plus_k) <- "numeric"
  
  # reconvert to matrix form
  bands_plus_k <- as.matrix(bands_plus_k)

  bands <- bands_plus_k[1:length(E),] # band cross section spectrum
  bands_lci <- bands_plus_k[(length(E)+1):(2*length(E)),] # band cross section spectrum CI
  bands_uci <- bands_plus_k[(2*length(E)+1):(3*length(E)),] # band cross section spectrum CI
  
  # generate full spectrum by adding up all bands contributions
  if(ncol(band_par)!=1){ # if there are multiple bands
    abs_cross <- rowSums(bands) # eq. 3 in publication
    abs_cross_lci <- sqrt(rowSums(bands_lci^2)) # eq. 5 in publication
    abs_cross_uci <- sqrt(rowSums(bands_uci^2)) # eq. 5 in publication
  } else {# if there is only one
    abs_cross <- bands
    abs_cross_lci <- bands_lci
    abs_cross_uci <- bands_uci
  }
  
  # optimal GMM model parameters
  k_opt <- bands_plus_k[(n_rows-1),]
  model_opt <- band_par[n_rows,]
  
  # return bands, band errors, spectrum, spectrum errors, and GMM model hyperparameters
  return(list(bands, bands_lci, bands_uci, abs_cross, abs_cross_lci, abs_cross_uci, model_opt, k_opt))
}

##########################################################################################################
###### functions to save absorption cross section full spectrum and transitions spectra
##########################################################################################################

save_sigma_eV_auto <- function(sigma, n_geoms, d, outlier_geoms, outlier_flag){
  # make directory to store results if doesn't exist
  outputDIR <- input_folder
  #if (!dir.exists(outputDIR)) {dir.create(outputDIR)}
  # set output file name
  file_name <- paste(outputDIR,'auto_d_spectra_eV_',n_geoms,'_geoms_', molecule,'.csv', sep = '')
  # set file header
  line <- paste('##### Reconstructed absorption cross section spectrum (',molecule,') ##### \n',
                '## Number of input transitions: ', n_states, '\n',
                '## Number of input geometries: ', n_geoms, '\n',
                '## Possible outliers (geometry id): ', paste(outlier_geoms, collapse = ', '), '\n',
                '## Remove outliers?: ', toString(outliers_flag), '\n',
                '## Model: auto-d \n',
                '## Empirical bandwidths d: ', paste(round(d,3), collapse =', '), '\n',
                '## Date: ', Sys.Date(), sep = '')
  
  # write header and dataframe to file
  write(line,file_name)
  suppressWarnings(
    write.table(
      format(sigma, digits = 6), file_name, sep = ",", row.names = FALSE, quote = FALSE, append = T
      )
    )
}

save_sigma_nm_auto <- function(sigma, n_geoms, d, outlier_geoms, outlier_flag){
  # eV to nm conversion
  c=299792458.0 # speed of light m/s
  hev=4.135667662e-15 # Planck's constant eV·s
  factor=hev*c*1.0e9
  
  sigma[,1] <- factor/sigma[,1]
  
  # change column name
  colnames(sigma)[1] <- 'wl_nm'
  
  # make directory to store results if doesn't exist
  outputDIR <- input_folder
  #if (!dir.exists(outputDIR)) {dir.create(outputDIR)}
  # set output file name
  file_name <- paste(outputDIR,'auto_d_spectra_nm_',n_geoms,'_geoms_', molecule,'.csv', sep = '')
  # set file header
  line <- paste('##### Reconstructed absorption cross section spectrum (',molecule,') ##### \n',
                '## Number of input transitions: ', n_states, '\n',
                '## Number of input geometries: ', n_geoms, '\n',
                '## Possible outliers (geometry id): ', paste(outlier_geoms, collapse = ', '), '\n',
                '## Remove outliers?: ', toString(outliers_flag), '\n',
                '## Model: auto-d \n',
                '## Empirical bandwidths d: ', paste(round(d,3), collapse =', '), '\n',
                '## Date: ', Sys.Date(), sep = '')
  
  # write header and dataframe to file
  write(line,file_name)
  suppressWarnings(
    write.table(
      format(sigma, digits = 6), file_name, sep = ",", row.names = FALSE, quote = FALSE, append = T
    )
  )
}

save_sigma_eV_gmm <- function(sigma, n_geoms, model_names, k_opt, outlier_geoms, outliers_flag){
  # make directory to store results if it doesn't exist
  outputDIR <- input_folder
  #if (!dir.exists(outputDIR)) {dir.create(outputDIR)}
  # set output file name
  file_name <- paste(outputDIR,'GMM_NEA_spectra_eV_',n_geoms,'_geoms_', molecule,'.csv', sep = '')
  # set file header
  line <- paste('##### Reconstructed absorption cross section spectrum (',molecule,') ##### \n',
                '## Number of input transitions: ', n_states, '\n',
                '## Number of input geometries: ', n_geoms, '\n',
                '## Possible outliers (geometry id): ', paste(outlier_geoms, collapse = ', '), '\n',
                '## Remove outliers?: ', toString(outliers_flag), '\n',
                '## Model: GMM-NEA \n',
                '## Number of mixtures (K): ', paste(k_opt, collapse = ', '), '\n',
                '## Model constraints (M): ', paste(model_names, collapse = ', '), '\n',
                '## Date: ', Sys.Date(), sep = '')
  
  # write header and dataframe to file
  write(line,file_name)
  suppressWarnings(
    write.table(
      format(sigma, digits = 6), file_name, sep = ",", row.names = FALSE, quote = FALSE, append = T
      )
    )
}

save_sigma_nm_gmm <- function(sigma, n_geoms, model_names, k_opt, outlier_geoms, outliers_flag){
  # eV to nm conversion
  c=299792458.0 # speed of light m/s
  hev=4.135667662e-15 # Planck's constant eV·s
  factor=hev*c*1.0e9
  
  sigma[,1] <- factor/sigma[,1]
  
  # change column name
  colnames(sigma)[1] <- 'wl_nm'
  
  # make directory to store results if it doesn't exist
  outputDIR <- input_folder
  #if (!dir.exists(outputDIR)) {dir.create(outputDIR)}
  # set output file name
  file_name <- paste(outputDIR,'GMM_NEA_spectra_nm_',n_geoms,'_geoms_', molecule,'.csv', sep = '')
  # set file header
  line <- paste('##### Reconstructed absorption cross section spectrum (',molecule,') ##### \n',
                '## Number of input transitions: ', n_states, '\n',
                '## Number of input geometries: ', n_geoms, '\n',
                '## Possible outliers (geometry id): ', paste(outlier_geoms, collapse = ', '), '\n',
                '## Remove outliers?: ', toString(outliers_flag), '\n',
                '## Model: GMM-NEA \n',
                '## Number of mixtures (K): ', paste(k_opt, collapse = ', '), '\n',
                '## Model constraints (M): ', paste(model_names, collapse = ', '), '\n',
                '## Date: ', Sys.Date(), sep = '')
  
  # write header and dataframe to file
  write(line,file_name)
  suppressWarnings(
    write.table(
      format(sigma, digits = 6), file_name, sep = ",", row.names = FALSE, quote = FALSE, append = T
    )
  )
}

##########################################################################################################
### function to plot spectra and save figures
##########################################################################################################

plot_spectra_eV <- function(sigma_auto_d,sigma_gmm,molecule){
  par(mar=c(4.5,4.6,1.5,1)+.1) # 'bottom', 'left', 'top', 'right'
  y_max <- 1.05*max(c(sigma_auto_d[,2],sigma_gmm[,2]))
  
  E <- sigma_auto_d[,1]
  
  plot(E,sigma_auto_d[,2], 
       xlab = 'E (eV)', 
       ylab= expression(sigma[abs] *'(E) (cm'^2*')'), 
       ylim = c(0,y_max), 
       type = "l", 
       lty = 1, lwd = 2,
       main = paste('NEA spectra for ', n_geoms, ' geometries', sep = ''),
       cex.axis = 1.2,
       cex.lab = 1.5)
  
  lines(E,sigma_gmm[,2], col = 2, lwd = 2)
  legend('top',legend = c(expression('auto-' ~ delta), 'GMM-NEA'), col = c(1,2), lty = c(1,1), lwd = 2, cex = 0.8)
  
  # Fill area between CI lines
  ul <- sigma_auto_d[,2]+sigma_auto_d[,4]  
  ll <- sigma_auto_d[,2]-sigma_auto_d[,3]
  polygon(c(E, rev(E)), c(ul, rev(ll)),
          col = scales::alpha(1, 0.1), lty = 0)  
  
  ul <- sigma_gmm[,2]+sigma_gmm[,4]  
  ll <- sigma_gmm[,2]-sigma_gmm[,3]
  polygon(c(E, rev(E)), c(ul, rev(ll)),
          col = scales::alpha(2, 0.1), lty = 0)
}

plot_spectra_nm <- function(sigma_auto_d,sigma_gmm,molecule){
  par(mar=c(4.5,4.6,1.5,1)+.1) # 'bottom', 'left', 'top', 'right'
  y_max <- 1.05*max(c(sigma_auto_d[,2],sigma_gmm[,2]))
  
  # eV to nm conversion factor
  c=299792458.0 # speed of light m/s
  hev=4.135667662e-15 # Planck's constant eV·s
  factor=hev*c*1.0e9
  
  wl <- factor/sigma_auto_d[,1]
  
  plot(wl,sigma_auto_d[,2], 
       xlab = expression(lambda *' (nm)'), 
       ylab= expression(sigma[abs] * '(' * lambda *') (cm'^2*')'), 
       ylim = c(0,y_max), 
       type = "l", 
       lty = 1, lwd = 2,
       main = paste('NEA spectra for ', n_geoms, ' geometries', sep = ''),
       cex.axis = 1.2,
       cex.lab = 1.5)
  
  lines(wl,sigma_gmm[,2], col = 2, lwd = 2)
  legend('top',legend = c(expression('auto-' ~ delta), 'GMM-NEA'), col = c(1,2), lty = c(1,1), lwd = 2, cex = 0.8)
  
  # Fill area between CI lines
  ul <- sigma_auto_d[,2]+sigma_auto_d[,4]  
  ll <- sigma_auto_d[,2]-sigma_auto_d[,3]
  polygon(c(wl, rev(wl)), c(ul, rev(ll)),
          col = scales::alpha(1, 0.1), lty = 0)  
  
  ul <- sigma_gmm[,2]+sigma_gmm[,4]  
  ll <- sigma_gmm[,2]-sigma_gmm[,3]
  polygon(c(wl, rev(wl)), c(ul, rev(ll)),
          col = scales::alpha(2, 0.1), lty = 0)
}


save_plots <- function(sigma_auto_d,sigma_gmm,molecule,n_geoms){
  # save plots in eV
  pdf(paste(input_folder,'Spectra_eV_auto_d_vs_GMM_NEA_',n_geoms,'_geoms_', molecule,'.pdf', sep = ''),
      width = 8.5, height = 6)
  plot_spectra_eV(sigma_auto_d, sigma_gmm, molecule)
  garbage <- dev.off()
  
  # save plots in nm
  pdf(paste(input_folder,'Spectra_nm_auto_d_vs_GMM_NEA_',n_geoms,'_geoms_', molecule,'.pdf', sep = ''),
      width = 8.5, height = 6)
  plot_spectra_nm(sigma_auto_d, sigma_gmm, molecule)
  garbage <- dev.off()
}