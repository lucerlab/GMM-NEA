#!/usr/bin/env Rscript

#######################################################################
#######################################################################
#####                                                             #####
#####     CODE TO GENERATE NEA SPECTRA USING AUTO-D AND NEA-GMM   #####
#####                                                             #####
#######################################################################
#######################################################################

#######################################################################
###########     STORE ARGUMENTS PASSED FROM TERMINAL     ##############
#######################################################################

# terminal example: Rscript GMM_NEA.R target_uracile.csv U6OH outliers ci
args <- commandArgs(trailingOnly=TRUE)

# test if all needed arguments are passed: if not, return an error
if (length(args)==1) {
  stop("At least two argument must be supplied (input file, output name).", call.=FALSE)
} 

# store input file name and molecule name
input_file <- args[1]
molecule <- args[2]

# by default, outliers are computed, confidence intervals are NOT
outliers <- TRUE
ci <- FALSE
if(any(args == 'no_outliers')){ outliers <- FALSE }
if(any(args == 'ci')){ ci <- TRUE }


######################################################################
###########    LOAD PACKAGES AND AUXILIARY FUNCTIONS    ##############
######################################################################

# set random seed
set.seed(1)

# set working directory
#setwd(getwd())

# check for missing packages
needed.packages <- c('stats','scales','robustbase','matrixStats','mclust','foreach','doParallel','parallel')
missing.pkg <- FALSE
pkgs <- c()
for(ii in 1:length(needed.packages)){
  if(length(find.package(needed.packages[ii], quiet = TRUE))==0) {
    pkgs <- c(pkgs,ii)
    missing.pkg <- TRUE
  }
}
if(missing.pkg){
  message(needed.packages[pkgs])
  stop("The above packages are missing. Please, install.", call.=FALSE)
}

# load specific packages
list.of.packages <- c('mclust','foreach')
for(package.i in list.of.packages){
  suppressWarnings(
    suppressMessages(
      library(
        package.i,
        character.only = TRUE
      )
    )
  )
}

#### load file with auxiliary functions
source('GMM_NEA_helper_funcs.R', chdir = TRUE)

######################################################################
###########    SET-UP MACHINE FOR PARALLEL COMPUTING    ##############
######################################################################

# set number of cores to use
n.cores <- parallel::detectCores() - 1

# create the cluster
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

# check cluster definition (optional)
print(my.cluster)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# make the cluster aware of external functions
parallel::clusterExport(my.cluster, c('weighted.IQR','auto_d','band_gen','band_gen_pdf','gauss','gmm_band'))


######################################################################
###########    LOAD DATASET AND FIND/REMOVE OUTLIERS    ##############
######################################################################

# Load full QM dataset
df_tr = read.csv(input_file, sep = ",", dec = ".", comment.char = "#", na.strings=c("","NA"))

# store number of states (transitions)
n_states <- ncol(df_tr)/2

# store number of geometries/configurations
n_geoms <- nrow(df_tr)

# Determine outliers with Mahalanobis distance and the False Discovery Rate (FDR)
if(outliers){
  outliers <- find_outliers(df_tr,0.001)
  df_tr <- outliers[[1]] # clean data frame
  outlier_geoms <- outliers[[2]] # list of outlier geometries
} else {
  outlier_geoms <- 'None'
}

######################################################################
###########                COMPUTE SPECTRA              ##############
######################################################################

# Arrays to store results (two columns per state plus columns for 
# energy E_ev and full spectrum with error)

# generate column names
cols <- c('E_ev', 'sigma_cm2_full','sigma_cm2_full_ci')

for(ii in 1:n_states){
  cols <- c(cols,
            paste('sigma_cm2_Band_',ii, sep = ''),
            paste('sigma_cm2_Band_',ii, '_ci', sep = ''))
}

# define spectral coverage
# e_min <- floor(10*min(df_tr[,1:n_states]))/10
# e_max <- ceiling(10*max(df_tr[,1:n_states]))/10
e_min <- floor(min(df_tr[,1:n_states]))
e_max <- ceiling(max(df_tr[,1:n_states]))
E <- seq(e_min, e_max, by = 0.01)

# create data frames to store results
sigma_auto_d <- as.data.frame(matrix(0, nrow = length(E), ncol = (2*n_states+3)))
colnames(sigma_auto_d)<- cols
sigma_gmm <- sigma_auto_d

message('Running auto-d...')

# compute bands and spectrum using auto-d
spectra_auto_d <- auto_d_spectra(df_tr, E, ci)
sigma_auto_d[,1] <- E
sigma_auto_d[,seq(4,(2*n_states+3), by = 2)] <- spectra_auto_d[[1]]
sigma_auto_d[,seq(5,(2*n_states+3), by = 2)] <- spectra_auto_d[[2]]
sigma_auto_d[,2] <- spectra_auto_d[[3]]
sigma_auto_d[,3] <- spectra_auto_d[[4]]
d_opt <- spectra_auto_d[[5]]

message('Running GMM-NEA...')

# compute bands and spectrum using GMM-NEA
spectra_gmm <- gmm_nea_spectra(df_tr,E,ci)
sigma_gmm[,1] <- E
sigma_gmm[,seq(4,(2*n_states+3), by = 2)] <- spectra_gmm[[1]]
sigma_gmm[,seq(5,(2*n_states+3), by = 2)] <- spectra_gmm[[2]]
sigma_gmm[,2] <- spectra_gmm[[3]]
sigma_gmm[,3] <- spectra_gmm[[4]]
gmm_model_opt <- spectra_gmm[[5]]
gmm_k_opt <- spectra_gmm[[6]]

# stop using parallel cluster
parallel::stopCluster(cl = my.cluster)

######################################################################
###########            SAVE AND PLOT RESULTS            ##############
######################################################################

message('Saving results...')

# save dataframes with spectra
save_sigma_eV_auto(sigma_auto_d,n_geoms,d_opt,outlier_geoms)
save_sigma_eV_gmm(sigma_gmm,n_geoms,gmm_model_opt,gmm_k_opt,outlier_geoms)

# save plot
png(paste(molecule,'/Spectra_eV_auto_d_vs_GMM_NEA_',n_geoms,'_geoms_', molecule,'.png', sep = ''),
    width = 8.5, height = 6,
    units = 'in',
    res = 300)
plot_spectra(sigma_auto_d,sigma_gmm,molecule)
garbage <- dev.off()


message('Done!!')

#renv::dependencies()
