################################################################################
#   MODULE FOR CLUSTERING FACTOR SCORES                                        #
################################################################################

################################################################################
### STEP 5, CLUSTER FACTOR SCORES
print('Running cluster.R')

t_cls0 <- proc.time()[3]                            # Mark time

t_cluster0 <- proc.time()[3] # Mark time

invisible(gc())

t_pop_bigfs1 <- 0                                  # Initialize timing stat
if (!exists('big.fs')) {
  big.fs <- attach.resource(paste0(opath_fs,       # Matrix of factor scores
		'fscore.desc')) }
if (!exists('idx_bigfs')){
  idx_bigfs <- mwhich(big.fs, 1, NA, 'neq')}		   # index of non NA rows

cat('Clustering...')
t_calc_cls0 <- proc.time()[3]                       # Mark time

cl <- start_par()                                  # Start parallelism

km <- bigkmeans(big.fs[idx_bigfs, 1:nfac], 		 # Cluster in parallel
		centers=nk,               
		iter.max = iter_trial,
		nstart = iter_seed,
		dist = 'euclid')

stop_par()                                         # Stops paral. if in Windows

t_calc_cls1 <- proc.time()[3] - t_calc_cls0          # Mark time
helap_calc_cls <- sprintf('%0.4f', t_calc_cls1/60^2) # Hrs elapsed
cat(paste('Done. Hrs elapsed:', helap_calc_cls, '\n'))

# Plot
#hist( rep(1:length(km$size), km$size))
#barplot(km$size, names.arg=1:length(km$size))

t_wr_rds0 <- proc.time()[3]                        # Mark time

cl.px <- rep(NA, nrow(big.fs))			# get cluster assgnmts
cl.px[idx_bigfs] <- km$cluster			# keeping NA structure from big.fs
saveRDS(cl.px, paste0(opath_cls,'cluster_pixels.rds'))

saveRDS(km$cluster, file=paste0(opath_cls,        	 # Save cluster number
				'cluster.rds')) 		   #  assignments
saveRDS(km$center, file=paste0(opath_cls,        	 # Save centriods
				'centers.rds'))
km[which(names(km) %in% 'cluster')] <- NULL        # Remove vec of cls nums
saveRDS(km, file=paste0(opath_cls,                 # Save without vector of
			'km_results.rds'))         #  cluster numbers
t_wr_rds1 <- proc.time()[3] -                      # Accumul. time
             t_wr_rds0 + t_wr_rds1

# save cluster assignments as factor centroid values.
dtc <- data.table(km$center)				# and by cluster
dtc[, cl := .I]
dtk <- setDT(list(cl.px))
kf <- merge(dtk, dtc, all=TRUE, by.x='V1', by.y='cl', sort=FALSE)
names <- paste0("f", 2:ncol(kf),".centroid")
names(kf) <- c("cl",names)
for(J in 2:ncol(kf)){
	file.name <- paste0(opath_cls, "f",J-1,"centroids.rds")
	saveRDS(kf[[J]], file=file.name)
	}

rm(big.fs,kf,dtc,dtk)
invisible(gc())
t_wr_ras0 <- proc.time()[3]                        # Mark time

# Write raster output
if (mk_ras == TRUE) {
	print('Writing raster output')
	l.files <- list.files(opath_cls, pattern = "\\.rds$", full.names=TRUE)
	l.files <- l.files[-grep('centers|results|cluster.rds', l.files)]    # Exclude centers and km_results files
		lapply(l.files,
		function(x) {annual.rasters(annual.data=x,
            	data.name='cls', yrs=npy, save.dir=opath, rtmpdir=ras.tmpdir, 
			remove.hrs=flushtemps.hrs, 
			template=paste0(srcpath,'mask.tif'))})
		}
t_wr_ras1 <- proc.time()[3] - t_wr_ras0 + t_wr_ras1# Accumul. time

rm(cl.px); invisible(gc())

restart <- data.frame(pr_complete = TRUE,          # Preproc complete (boolean)
           fs_complete = TRUE,                     # FacAnal complete (boolean)
           cl_complete = TRUE,                     # Cluster complete (boolean)
           I = NA, J = NA)                         # Copies of loop vals
t_wr_rds0 <- proc.time()[3]                        # Mark time
saveRDS(restart, paste0(bpath,'restart.rds'))      # Save restart values
t_wr_rds1 <- proc.time()[3] -                      # Accumul. time
             t_wr_rds0 + t_wr_rds1
t_cls1 <- t_cls0 - proc.time()[3]                    # Mark time
helap_cls <- sprintf('%0.4f', t_fs1/60^2)           # Hrs elapsed
print('cluster.R Done')
invisible(gc())
