################################################################################
#   GLOBAL VARIABLES                                                           #
################################################################################

cat('Loading libraries & variables...')
### Load libraries
os_env <- Sys.info()[['sysname']]                  # Query OS type

libs <- c('raster', 'parallel', 'foreach', 'biganalytics', 'zoo',
	   'PolarMetrics', 'ggplot2', 'data.table', 'expm')
# Add libs for shared memory based on OS type
if (os_env == 'Linux') {
  libs <- c(libs, 'doMC') 
} else if (os_env == 'Windows') {
  libs <- c(libs, 'doSNOW', 'snow', 'doParallel')
}
checkLoad(libs)                                    # Check for, install & load
options(bigmemory.typecast.warning=FALSE)          # Stop typecast warnings

ncores <- detectCores() - 1                        # Set no. cores to use
Zth <- 0.75                                        # Pppn of pix-yr allowed to
                                              #  contain zeros & pass imputing
Zth2 <- 0.25                                       # Pppn of total pix vals
                                      #  allowed contain zeros & pass imputing

nk=50                                              # Num of desired clusters
### Diagnostics variables (queried for HTML page)
arch <- Sys.info()[['machine']]                    # Query mach. architecture
user <- Sys.info()[['user']]                       # Query name
rtme <- format(Sys.time(), "%Y-%m-%d-%H:%M")       # Time
run_id <- paste(rtme, user, os_env, arch,          # Run identifier
                paste0(ncores, 'cores'),
                paste0(nk, 'clusters'), sep='_')

### Directory paths
bpath <- 'F:/'                                     # Base path
hpath <- paste0(hpath, 'landat/'                   # Home path (this project)
#hpath <- '/run/media/bjorn/phenometric/LanDAT_code_share/' # Base path
opath <- paste0(hpath, 'output/')                  # Output data dir
opath_pm <- paste0(opath, 'pmetric/')              # Output dir for polar met.
opath_fs <- paste0(opath, 'fscore/')               # Output dir for factor sc.
opath_cls <- paste0(opath, 'cluster/')             # Output dir for clust asgn.
opath_it <- paste0(opath, 'it/')                   # Output dir IT data
opath_ras <- paste0(opath, 'raster_results/')      # Output dir for rasters
srcpath <- paste0(hpath, 'data/tx_subset/')        # Original NDVI source
dgpath <- paste0(opath, 'diagnostics/')            # Path to diagnostics dir


### Clean up or create output directories
if (dir.exists(opath)) {
  run_id_back <- readRDS(paste0(opath, 'run_id.rds'))
  run_id_back <- sub(':', '', run_id_back)
  if (os_env == 'Windows') {
    files_fr <- list.files(opath, full.names=TRUE)
    files_to <- sub('output\\/',
                    paste0('output_back/', run_id_back, '/'), files_fr)
    dir.create(sub('output\\/',     # Create output dir
                    paste0('output_back/', run_id_back, '/'), opath),
               recursive=TRUE)
    msg1 <- file.rename(from=files_fr,             # Move (not copy) files
                 to=files_to)                      #  by renaming them
  } else {
    msg1 <- file.rename(opath, paste0(hpath,       # Rename prior output dir
                        'output_', run_id_back))
  }
  if (msg1[1] == FALSE) {
    stop('vars.R (msg1): failed to rename opath directory')
  }
  rm(msg1)
}
dir.create(opath_pm, recursive=TRUE)               # Create output dirs
dir.create(opath_fs, recursive=TRUE)
dir.create(opath_cls, recursive=TRUE)
dir.create(opath_it, recursive=TRUE)
msg2 <- file.copy(from=paste0(hpath,               # Copy diagnostics source
			    'diagnostics/'),       #  dir to output dir
          to=opath, recursive=TRUE)
if (msg2 == FALSE) {
  stop('vars.R (msg2): failed to copy diagnostics/')
}
rm(msg2)
saveRDS(run_id, file=paste0(opath, 'run_id.rds')) # Label output

### Global otions
continue <- FALSE                                  # Cont. from prior run end?
mask <- TRUE                                       # Use mask raster?
mk_ras <- TRUE                                     # Save intermed. output as
                                                   #  rasters

### Preprocessing options
maxcells <- 10^4                                   # Pix / chunk (preprocess.R)
pptrn <- 'A20[0-9][0-9].tif$'                        # Input NDVI file name ptrn
f2pr <- list.files(srcpath,                        # Rasters to ingest
		   pattern=pptrn, full.names=TRUE)
if (length(f2pr) == 0) {
  stop('vars.R (f2pr): Script stopped. No rasters to process')
}
# Poll raster mask for dimension variables
ras <- brick(f2pr[1])
if (mask == TRUE) { # Use prexisting mask to limit num pixels to process
  ras_mask <- raster(paste0(srcpath, 'mask.tif'))  # Set raster file for mask
  cid2pr <- Which(ras_mask, cells=TRUE)            # Raster cellIDs to procs.
  ncell <- length(cid2pr)
} else {
  ncell <- ncell(ras)                              # Num pix to process
}

saveRDS(ncell, file=paste0(opath, 'ncell.rds'))    # Save num pixels
nch <- ceiling(ncell/maxcells)                     # Num chunks loop over
tmp <- intv_starts(cid2pr, n=nch)                  # Find breakpoints & freq
cid_starts <- tmp[[1]]                             # Cellid start points
cid_freq <- tmp[[2]]                               # Num cells in each chunk
rm(tmp)

spc <- nlayers(ras)                                # Tot num samples per pix
nsamp <- spc * length(f2pr)                        # Tot num samples per pix
years <- as.integer(sub('\\.tif$', '',             # Years of NDVI data
			sub('.*\\/y2', '2', f2pr)))
nyr <- length(years)                               # Tot num yrs data in input
npy <- nyr-1                                       # Num phen yrs (PhenoMetr)
yrbeg <- min(years)                                # first yr of NDVI data
yrend <- max(years)                                # last yr of NDVI data
doy <- rep(seq(from=3, to=363, by=8), nyr)         # Time vals for NDVI input
pfx <- paste(yrbeg, yrend, sep='-')                # NDVI input file prefix


### Factor analysis options
lb_unique <- 0.05                                  # FctAnl lw bnd for uniqnes
nfac <- 4                                          # Num of factor dimensions

### Clustering options
iter_seed <- 10^0                                  # Itrs to find cl seed vals
iter_trial <- 10^0                                 # Max no. seed trials
minobsyr <- 2                                      # Min num complete yrs
                                                   #  input needed to GF

### IT options
window <- 5					   # Mov.win.dim. along 1side
year.range <- 1:npy

### Raster processing options
rast.format <- 'GTiff'				   # File format, e.g., 'GTiff'
ras.tmpdir <- paste0(bpath, 'tmp'		   # Tmp. dir for ras process
flushtemps.hrs <- 0.3				   # Num hrs before flushing raster temporary dir (see raster::removeTmpFiles)

### Miscillaneous variables
ram_avail_kb <- checkRAM(os_env)                   # Check available RAM (kB)
# Poll calc_metrics() variable names
pmvarnames <- colnames(calc_metrics(1:4, t=rep(1:2,times=2),
						  yr_type='cal_yr', spc=2,
						  lcut=0.15, hcut=0.8,
						  return_vecs=F, sin_cos=T))

if (continue == TRUE) {
  restart <- readRDS(paste0(hpath,'restart.rds'))  # Load restart values
  Ibeg <- restart$I                                # Init Ibeg to rstrt value
  Jbeg <- restart$J                                # Init Jbeg to rstrt value
} else {
  restart <- data.frame(pr_complete = FALSE,       # Preproc complete (boolean)
             fs_complete = FALSE,                  # FacAnal complete (boolean)
             cl_complete = FALSE,                  # Cluster complete (boolean)
             I = NA, J = NA)                       # Copies of loop vals
  Ibeg <- 1                                        # Init Ibeg to 1
  Jbeg <- 1                                        # Init Jbeg to 1
}
### Example program restart codes for interruption during fscores.R
### pr_complete = TRUE,
### fs_complete = FALSE,
### cl_complete = FALSE,
### I = 1, J = NA
### Means prog has has begun 1st I loop, but has not begun J loop yet.
### Note that parallel loops (P1 or P2) can not be restarted and are
### not included in this restart procedure. The entire parallel loop must
### be redone.

# Initialize variables for diag.R
t_wr_rds1 <- 0                                     # Initialize timing stat
t_wr_ras1 <- 0                                     # Initialize timing stat
helap_impute <- 0
helap_calc_pm <- 0
helap_calc_fs <- 0
helap_calc_cls <- 0
helap_calc_it <- 0
helap_wr_rds_pm <- 0
helap_wr_ras_pm <- 0
helap_pop_bigpm <- 0
helap_calc_pm <- 0

cat('Done.\n') 
