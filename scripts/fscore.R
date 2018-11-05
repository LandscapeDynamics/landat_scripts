################################################################################
#   MODULE FOR CALCULATING FACTOR SCORES                                       #
################################################################################

################################################################################
### STEP 4, CALCULATE FACTOR SCORES
print('Running fscore.R')

t_fs0 <- proc.time()[3]                            # Mark time

t_load_pm1 <- 0                                    # Initialize timing stat
t_pop_bigpm1 <- 0                                  # Initialize timing stat

invisible(gc())
if (!exists('big.pm')) {
  big.pm <- attach.big.matrix(paste0(opath_pm,     # Attach to NDVI data
				     'pmetric.desc'))
}

# index of rows in big.pm with no NAs
idx_bigpm <- mwhich(big.pm, 1:ncol(big.pm), NA, 'neq', 'OR')	# QUESTION: are NA pixels alwasy NA for all years?
ncellOK <- length(idx_bigpm)/npy						# If so, this will work. If not, need different method.	

saveRDS(idx_bigpm, paste0(opath,'idx_bigpm.rds'))			# save good pixel IDs

big.fs <- filebacked.big.matrix(nrow=nrow(big.pm),  # Init. file
				ncol=nfac,
				dimnames=list(rownames=NULL,
					      colnames=paste0('factor', 1:nfac)),
				backingfile='fscore.bigmatrix',
				backingpath=opath_fs,
				descriptorfile='fscore.desc',
				init=NA)


# Scale each incoming polar metric variable first before factor analysis
#  to accelerate computation
for (Z in 1:ncol(big.pm)) {
  big.pm[idx_bigpm, Z] <- normalize_center_sint(big.pm[idx_bigpm,
						Z])# Center & normalize to fill
}
                                                   #  range of data type short

# Factor Analysis of pheno-metrics
# Use the equation form in factanal b/c it allows for explicit control of NA's
#  eqn <- paste('~', paste(vnames[-1], collapse='+'), sep='')
#  eqn <- eval(parse(text=paste('~x$', paste(vnames[-1], collapse='[]+x$'), '[]', sep='')))
#  fit <- factanal(eqn, 4, rotation='varimax', scores='Bartlett',
#                  lower=0.005, na.action=na.exclude) # Only process non-NA
#NC <- ncol(pm) # No. of columns

cat('Running factor analysis...')
t_calc_fa0 <- proc.time()[3]                       # Mark time
### Run factor analysis
fa <- factanal(big.pm[idx_bigpm, 1:ncol(big.pm)], # Exploratory FA
	       nfac, rotation='varimax', lower=lb_unique)
t_calc_fa1 <- proc.time()[3] - t_calc_fa0          # Num seconds elapsed
helap_calc_fa <- sprintf('%0.4f', t_calc_fa1/60^2) # Hrs elapsed
cat(paste('Done. Hrs elapsed:', helap_calc_fa, '\n'))
t_wr_rds0 <- proc.time()[3]                        # Mark time
saveRDS(fa, file=paste0(opath_fs, 'fa_results.rds'))
t_wr_rds1 <- proc.time()[3] -                      # Accumul. time
             t_wr_rds0 + t_wr_rds1


pm.sds <- colsd(big.pm, na.rm=TRUE)

big.fs[idx_bigpm,1:nfac] <- t(t(big.pm[idx_bigpm, 1:ncol(big.pm)])/pm.sds) %*% # Regr. factor
            solve(fa$correlation) %*% loadings(fa) #  score (Thomson's method)

idx_bigfs <- mwhich(big.fs, 1, NA, 'neq')		   # index of non NA rows
saveRDS(idx_bigfs, paste0(opath_fs, "bigfs_goodpixelIDs.rds"))

invisible(gc())


t_wr_ras0 <- proc.time()[3]                        # Mark time
# Write raster output
if (mk_ras == TRUE) {
  print('Writing raster output')
  factnames <- paste0("fscore", 1:nfac)
  t_wr_rds0 <- proc.time()[3]                      # Mark time

    annual.rasters.m(annual.data=big.fs,
    data.name='fscore', var.names=factnames,
    yrs=npy, template=paste0(srcpath, 'mask.tif'),
    datatype='FLT4S', save.dir=opath, rtmpdir=ras.tmpdir, remove.hrs=flushtemps.hrs)
  t_wr_rds1 <- proc.time()[3] -                    # Accumul. time
               t_wr_rds0 + t_wr_rds1
}
t_wr_ras1 <- proc.time()[3] - t_wr_ras0 + t_wr_ras1# Accumul. time

restart <- data.frame(pr_complete = TRUE,          # Preproc complete (boolean)
           fs_complete = TRUE,                     # FacAnal complete (boolean)
           cl_complete = FALSE,                    # Cluster complete (boolean)
           I = NA, J = NA)                         # Copies of loop vals
t_wr_rds0 <- proc.time()[3]                        # Mark time
saveRDS(restart, paste0(bpath,'restart.rds'))      # Save restart values
t_wr_rds1 <- proc.time()[3] -                      # Accumul. time
             t_wr_rds0 + t_wr_rds1
t_fs1 <- t_fs0 - proc.time()[3]                    # Mark time
helap_fs <- sprintf('%0.4f', t_fs1/60^2)           # Hrs elapsed
print('fscore.R Done')
invisible(gc())
