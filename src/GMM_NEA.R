#!/soft/R-4.1.1/bin/Rscript

# Author: Luis Cerdán <lcerdanphd_at_gmail.com>
#                     <https://github.com/lucerlab>

#######################################################################
#######################################################################
#####                                                             #####
#####     CODE TO GENERATE NEA SPECTRA USING AUTO-D AND NEA-GMM   #####
#####                                                             #####
#######################################################################
#######################################################################

######################################################################
###########    LOAD PACKAGES AND AUXILIARY FUNCTIONS    ##############
######################################################################

# set random seed
set.seed(1)

# check for missing packages
needed.packages <- c('stats','scales','robustbase','matrixStats','mclust','foreach','doParallel','parallel','optparse')
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

# load specific packages. Rest of packages will be called using :: notation
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


#######################################################################
###########     STORE/ASSESS ARGUMENTS PASSED FROM TERMINAL     #######
#######################################################################

#library("optparse")

# terminal example:
# Rscript GMM_NEA.R -p "uracil radical" -f vee_f_uracil_radical.csv -m U6OH [-o TRUE] [-c FALSE] [-n -1]
# or
# Rscript GMM_NEA.R --path "uracil radical" -f vee_f_uracil_radical.csv -m U6OH [-o TRUE] [-c FALSE] [-n -1]

option_list = list(
  optparse::make_option(c("-p", "--path"), type="character", default=NULL, 
              help="Path to folder with dataset [default=current folder]", metavar="character"),
  optparse::make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Dataset file name", metavar="character"),
  optparse::make_option(c("-m", "--molecule"), type="character", default="unknown", 
              help="Output file name identifier [default=%default]", metavar="character"),
  optparse::make_option(c("-t", "--threads"), type="integer", default= -1, 
              help="Number of threads [default=%default]", metavar="integer"),
  optparse::make_option(c("-o", "--outliers"), type="logical", default= TRUE, 
              help="Find and remove outliers? [default=%default]", metavar="logical"),
  optparse::make_option(c("-c", "--conf_inter"), type="logical", default= FALSE, 
              help="Compute confidence intervals? [default=%default]", metavar="logical"),
  optparse::make_option(c("-l", "--lowerE"), type="double", default= NULL, 
              help="Lower energy bound for spectra (eV) [default=%default]", metavar="double"),
  optparse::make_option(c("-u", "--upperE"), type="double", default= NULL, 
              help="Upper energy bound for spectra (eV) [default=%default]", metavar="double"),
  optparse::make_option(c("-d", "--deltaE"), type="double", default= 0.01, 
              help="Energy step size for spectra (eV) [default=%default]", metavar="double")
); 

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);


# test if all needed arguments are passed: if not, return an error
if ( is.null(opt$file) ) {
  stop("The argument '--file' (-f) must be supplied.", call.=FALSE)
}

# set new variables with the arguments
if ( is.null(opt$path) ) { # if path is NULL, set to current folder
  input_folder <- ''
} else {
  input_folder <- paste(opt$path, '/', sep = '')
}
input_file <- opt$file
molecule <- opt$molecule
outliers <- opt$outliers
ci <- opt$conf_inter
n_threads <- opt$threads


######################################################################
###########    SET-UP MACHINE FOR PARALLEL COMPUTING    ##############
######################################################################

# get number of available threads
n.cores <- parallel::detectCores()

# If n_threads = -1, use all, else, use minimum of n_threads and n.cores
if ( n_threads != -1) { n.cores <- min(c(n.cores,n_threads))}

# create the cluster
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK" # if working in Linux, you can switch "PSOCK" to "FORK"
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

# Load full QM dataset (VEE and f)
df_tr = read.csv(paste(input_folder, input_file, sep = ''), sep = ",", dec = ".", comment.char = "#", na.strings=c("","NA"))

# store number of states/transitions
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

# Arrays to store results (three columns per state plus columns for 
# energy E_ev and full spectrum with errors)

# generate column names
cols <- c('E_ev', 'sigma_cm2_full','sigma_cm2_full_lci','sigma_cm2_full_uci')

for(ii in 1:n_states){
  cols <- c(cols,
            paste('sigma_cm2_Band_',ii, sep = ''),
            paste('sigma_cm2_Band_',ii, '_lci', sep = ''),
            paste('sigma_cm2_Band_',ii, '_uci', sep = ''))
}

# define spectral coverage (check arguments pass to terminal)
if ( is.null(opt$lowerE) ) {
  e_min <- floor(min(df_tr[,1:n_states]))
} else {
  e_min <- opt$lowerE
}
if ( is.null(opt$upperE) ) {
  e_max <- ceiling(max(df_tr[,1:n_states]))
} else {
  e_max <- opt$upperE
}

dE_eV <- opt$deltaE

E <- seq(e_min, e_max, by = dE_eV)

# create data frames to store results
sigma_auto_d <- as.data.frame(matrix(0, nrow = length(E), ncol = (3*n_states+4)))
colnames(sigma_auto_d)<- cols
sigma_gmm <- sigma_auto_d

message('Running auto-d...')

# compute bands and spectrum using auto-d
spectra_auto_d <- auto_d_spectra(df_tr, E, ci)
sigma_auto_d[,seq(5,(3*n_states+4), by = 3)] <- spectra_auto_d[[1]] # add bands
sigma_auto_d[,seq(6,(3*n_states+4), by = 3)] <- spectra_auto_d[[2]] # add bands lci
sigma_auto_d[,seq(7,(3*n_states+4), by = 3)] <- spectra_auto_d[[3]] # add bands uci
sigma_auto_d[,2] <- spectra_auto_d[[4]] # add spectrum
sigma_auto_d[,3] <- spectra_auto_d[[5]] # add spec lci
sigma_auto_d[,4] <- spectra_auto_d[[6]] # add spec uci
d_opt <- spectra_auto_d[[7]] # add optimized linewidths


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

sigma_gmm[,seq(5,(3*n_states+4), by = 3)] <- spectra_gmm[[1]] # add bands
sigma_gmm[,seq(6,(3*n_states+4), by = 3)] <- spectra_gmm[[2]] # add bands lci
sigma_gmm[,seq(7,(3*n_states+4), by = 3)] <- spectra_gmm[[3]] # add bands uci
sigma_gmm[,2] <- spectra_gmm[[4]] # add spectrum
sigma_gmm[,3] <- spectra_gmm[[5]] # add spec lci
sigma_gmm[,4] <- spectra_gmm[[6]] # add spec uci
gmm_model_opt <- spectra_gmm[[7]] # add optimized model names
gmm_k_opt <- spectra_gmm[[8]] # add optimized number of mixtures

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
pdf(paste(input_folder,'Spectra_eV_auto_d_vs_GMM_NEA_',n_geoms,'_geoms_', molecule,'.pdf', sep = ''),
    width = 8.5, height = 6)
plot_spectra(sigma_auto_d,sigma_gmm,molecule)
garbage <- dev.off()


message('Done!!')
