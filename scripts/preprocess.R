################################################################################
#   MODULE FOR IMPUTING MISSING NDVI AND CALCULATING POLAR METRICS             #
################################################################################

t_impute1 <- 0                                     # Initialize timing stat
t_calc_pm1 <- 0                                    # Initialize timing stat

print('Running preprocess.R')

################################################################################
### STEP 1, CREATE BIG DATA MATRICES
cat('Imputing...')

t_impute0 <- proc.time()[3]                        # Mark time


big.pm <- filebacked.big.matrix(nrow=ncell*npy,    # Initialize file
				ncol=length(pmvarnames)-length(c(1:6,13)),
				dimnames=list(rownames=NULL,
					      colnames=pmvarnames[-c(1:6,13)]),
				backingfile='pmetric.bigmatrix',
				backingpath=opath_pm,
				descriptorfile='pmetric.desc',
				type='short', init=NA)

big.pm_doy <- filebacked.big.matrix(nrow=ncell*npy, # Initialize file
				ncol=length(2:6),
				dimnames=list(rownames=NULL,
					      colnames=pmvarnames[2:6]),
				backingfile='pmetric_doy.bigmatrix',
				backingpath=opath_pm,
				descriptorfile='pmetric_doy.desc',
				type='short', init=NA)

big.pm_aavg <- filebacked.big.matrix(nrow=ncell*npy, # Initialize file
				ncol=1,
				dimnames=list(rownames=NULL,
					      colnames=pmvarnames[13]),
				backingfile='pmetric_aavg.bigmatrix',
				backingpath=opath_pm,
				descriptorfile='pmetric_aavg.desc',
				type='double', init=NA)

big.nImpVals <- filebacked.big.matrix(nrow=nch,     # Initialize file
				ncol=1,
				dimnames=list(rownames=NULL,
					      colnames='num_imputed_vals'),
				backingfile='nImpVals.bigmatrix',
				backingpath=opath_pm,
				descriptorfile='nImpVals.desc',
				type='double', init=NA)

################################################################################
### STEP 2, LOAD NDVI AND IMPUTE MISSING NDVI VALUES

cl <- start_par()                                  # Start parallelism
# I loop is to process PM
cid_good <- foreach (I = 1:nch, .export=c('mwhich'),
			      .packages=c('PolarMetrics','raster','bigmemory'),
			      .combine=c) %dopar% {

  # Within foreach loop, must (re-)attach to big.matrix to access data
  big.pm <- attach.big.matrix(paste0(opath_pm,     # Attach to NDVI data
				     'pmetric.desc'))
  big.pm_doy <- attach.big.matrix(paste0(opath_pm, # Attach to NDVI data
				     'pmetric_doy.desc'))
  big.pm_aavg <- attach.big.matrix(paste0(opath_pm, # Attach to NDVI data
				     'pmetric_aavg.desc'))
  big.nImpVals <- attach.big.matrix(paste0(opath_pm, # Attach to NDVI data
				     'nImpVals.desc'))

  if (I < nch) {
    cid2pr_s <- cid2pr[which(cid2pr >= cid_starts[I] # Ras cell IDs to process
			     & cid2pr < cid_starts[I+1])]
  } else { # If last chunk then get CIDs to end of raster
    cid2pr_s <- cid2pr[which(cid2pr >= cid_starts[I])] # Ras cell IDs to process
  }
  ncell2 <- length(cid2pr_s)                              # Number of pix to process
  ndvi <- matrix(nrow = ncell2, ncol = nsamp)        # Initialize matrix
  spc_seq <- seq(1, ncol(ndvi), by=spc)            # Seq. of yr start points

  for(J in 1:length(f2pr)){
    x <- brick(f2pr[J])                            # Make raster brick
    dat <- extract(x, cid2pr_s)                 # Extract NDVI values
    frstCol <- spc_seq[J]                          # Idx of frst col to extrct
    lastCol <- frstCol + spc - 1                   # Idx of last col to extrct
    ndvi[ , frstCol:lastCol] <- dat                # Populate matrix
  }

  ### Impute
  # Use impute_median() to fill missing values. Do not use a parallel loop
  #  over the top of this function, as it contains parallel loops
  #  within. Output is returned as a list with 3 variables.
  t_impute0 <- proc.time()[3]                      # Mark time
  imOut <- impute_median(ndvi, spc=spc, Zth=Zth) # Impute

  ndvi_imputed <- imOut$imputed                    # NDVI vals after imputing 
  big.nImpVals[I,] <- imOut$nImpVals                # Save num imputed elemnts

  # Note the objective of the next two processes (Zcpgf & NAcpgf) is to:
  #  1) Replace rows with "too many" zeros with NA
  #  2) Replace rows that still have NA's after imputing with NA
  #  Note that a row in ndvi_imputed icludes all values for a pixel

  Zcpgf <- which((rowSums(                 # Rows exceeding Zth2 ppn of 0s
                          ndvi_imputed == 0) >= nsamp*Zth2) == T)

  if (any(length(Zcpgf))) {
    ndvi_imputed[Zcpgf, ] <- NA                    # Set rows to NA
  }                                                #  zeros after g.f.

  # If any part of row still has NA after gap-filling, set entire row to NA
  NAcpgf <- mwhich(ndvi_imputed, 1:nsamp, NA, 'eq', 'OR') # Rows with NA's

  if (any(length(NAcpgf))) {
    ndvi_imputed[NAcpgf, ] <- NA                   # Set to NA
    pix2pr <- (1:ncell2)[-NAcpgf]                    # Rows (also pixels) to analyze
    cidOK <- cid2pr_s[-NAcpgf]                     # CellIDs that passed QC
  } else {
    pix2pr <- 1:ncell2
    cidOK <- cid2pr_s
  }
 
  t_impute1 <- proc.time()[3] - t_impute0 +        # Accumulate time
               t_impute1


################################################################################
### STEP 3, CALCULATE PHENO-VARIABLES (POLAR METRICS) FROM NDVI

  t_calc_pm0 <- proc.time()[3]                     # Mark time

  # Loop over each row and recombine rows of output
  if (I == 1) {
    a0 <- 0
  } else {
    a0 <- sum(cid_freq[1:(I-1)]) * npy     # big.matrix row offset
  }                                                #  for Ith iteration
  for (J in 1:length(pix2pr)) {                              # Loop over cells in Ith sbst

    a <- a0 + (pix2pr[J]-1) * npy + 1                      # beg row in big.m to popul.
    b <- a + npy - 1                               # ending row

    # NOTE: That this preserves the correct order in big.matrices that corresponds
    #  to the raster mask by skipping NA's (i.e., rows that could not be imputed.)

    if (J %in% pix2pr) { # If this pixel passed QC then calc PolarMetrics

      # Note that if the output of pm happens to contain more than npy rows of
      #  data those rows are truncated off.
      pm <- calc_metrics(ndvi_imputed[pix2pr[J], ],
                         t=doy, yr_type='cal_yr', spc=spc,
                         lcut=0.15, hcut=0.8,
                         return_vecs=FALSE, sin_cos=TRUE)[1:npy,]

      t_calc_pm1 <- proc.time()[3] - t_calc_pm0      # Accumul. time
                    + t_calc_pm1

      big.pm_aavg[a:b, ] <- as.matrix(pm[ , 13])     # Populate ann.avg NDVI vars
      big.pm_doy[a:b, ] <- as.matrix(pm[ , 2:6])     # Populate DOY timing vars

      # Convert sin/cos & integer NDVI metric (x100) to signed short integers
      pm[ , 8:12] <- apply(pm[ , 8:12], MARGIN=2,    # Scale mag's to reduce sz
                           FUN=function(x){as.integer(x * 10^2)})
      pm[ ,14:23] <- apply(pm[ ,14:23], MARGIN=2,    # Scale sin/cos to red. sz
                           FUN=function(x){as.integer(x * 10^4)})

      big.pm[a:b, ] <- as.matrix(pm[, -c(1:6,13)])   # vars for fact. anal.
    }
  }

  if (I %% ncores == 0) {
    file.remove(list.files(path=opath,             # Remove prior dummy files
			   pattern='%complete', full.names=TRUE))
  }
  cat('', file=paste0(opath,                       # Write dummy file to show
		      sprintf('%08.4f', 100 * I/nch), '%complete')) # progress

  cidOK                                            # Repeat last foreach concat.

}                                                  # End I loop
stop_par()                                         # Stops par if in Windws

saveRDS(cid_good, paste0(opath,'cid_good.rds'))		# save cell IDs of good pixels


cat('Done.\n')
invisible(gc())

t_wr_ras0 <- proc.time()[3]                        # Mark time
# Write raster output
if (mk_ras == TRUE) {
  print('Writing PolarMetrics raster output')
  t_wr_rds0 <- proc.time()[3]                      # Mark time
  annual.rasters.m(annual.data=big.pm,
      data.name='pmetric', var.names=colnames(big.pm),
      yrs=npy, template=paste0(srcpath, 'mask.tif'),
      datatype='INT2S', save.dir=opath)

  print('Writing PolarMetrics DOY raster output')
  annual.rasters.m(annual.data=big.pm_doy,
      data.name='pmetric', var.names=colnames(big.pm_doy),
      yrs=npy, template=paste0(srcpath, 'mask.tif'),
      datatype='INT2S', save.dir=opath)

  print('Writing Annual Average NDVI raster output')
  annual.rasters.m(annual.data=big.pm_aavg,
      data.name='pmetric', var.names=colnames(big.pm_aavg),
      yrs=npy, template=paste0(srcpath, 'mask.tif'),
      datatype='FLT4S', save.dir=opath)

  t_wr_rds1 <- proc.time()[3] -                    # Accumul. time
               t_wr_rds0 + t_wr_rds1
}
t_wr_ras1 <- proc.time()[3] - t_wr_ras0 + t_wr_ras1# Accumul. time

helap_impute <- sprintf('%0.4f', t_impute1/60^2)   # Hrs elapsed
helap_calc_pm <- sprintf('%0.4f', t_calc_pm1/60^2) # Hrs elapsed

restart <- data.frame(pr_complete = TRUE,          # Preproc complete (T/F)
           fs_complete = FALSE,                    # FacAnal complete (T/F)
           cl_complete = FALSE,                    # Cluster complete (T/F)
           I = NA, J = NA)                         # Copies of loop vals

t_wr_rds0 <- proc.time()[3]                        # Mark time
saveRDS(restart, paste0(bpath,'restart.rds'))      # Save restart values
t_wr_rds1 <- proc.time()[3] -                      # Accumul. time
             t_wr_rds0 + t_wr_rds1

print('preproces.R Done')
invisible(gc())
