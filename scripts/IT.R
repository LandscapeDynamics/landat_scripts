################################################################################
#   MODULE FOR CALCULATING IT METRICS FROM CLUSTER OUTPUT                      #
################################################################################
#
# run full IT metric calculations
# many parameters required here are defined in the 'vars' script.
# IT functions required here are in the 'fxns' script.

### Use cluster.mean to get mean ndvi values by phenoclass.
	classes.path <- paste0(opath_cls, "cluster.rds")
	if(!exists('big.pm_aavg')) {
	big.pm_aavg <- attach.big.matrix(paste0(opath_pm, 	# attach to ndvi data
				     'pmetric_aavg.desc'))
	}
	idx_aavg <- mwhich(big.pm_aavg, 1, NA, 'neq')		# index of non-NA pixels = idx_bigpm
	metric.val <- big.pm_aavg[idx_aavg,1]			# a_avg variable
	pc.ndvi <- cluster.mean(pclasses=classes.path, metric.vals=metric.val)
	# save result
	saveRDS(pc.ndvi, paste0(opath_it, "clus_annualndvi_vals.rds"))

### load phenoclass attribute data
	factors <- data.table(readRDS(paste0(opath_cls,"centers.rds")))
		names(factors) <- paste0("f",1:ncol(factors))
	centroid.vals <- cbind(factors,pc.ndvi)		# used in IT.do

### use landscape.mean to get mean annual ndvi in landscape windows
	# set up paths to data
	rast.ndvi <- list.files(paste0(opath_ras,"pmetric/"), pattern="\\.tif$", 
		full.names=TRUE)
	rast.ndvi <- rast.ndvi[grep("a_avg", rast.ndvi)]

	land.ndvi <- landscape.mean(metric.vals=rast.ndvi, winsize=window, 	# use landscape.mean
		zRange=year.range, n.rows=25)	 
	# save result
	writeRaster(land.ndvi, filename=paste0(opath_it,"land",window,"_ndvi"), 
		format=rast.format, overwrite=T)
	meanndvi.path <- paste0(opath_it,"land",window,"_ndvi.tif")		# used in add.vars

### composite phenoclass year layers for data extraction.
	rast.cl <- list.files(paste0(opath_ras,"cls/"), 
		pattern="\\.tif$", full.names=TRUE)
	rast.cl <- rast.cl[grep("cluster", rast.cl)]
	cl.stack <- stack(rast.cl)
	writeRaster(cl.stack, filename=paste0(opath_it,"phenoclasses"), 
		format=rast.format, datatype='INT2U', overwrite=T)
	classes.path <- paste0(opath_it,"phenoclasses.tif")		# used in several IT functions

##### Execute IT, full run #####

# use IT.parRun to execute both IT.priors and IT.do on data
	IT.parRun(rasterData=classes.path, numClasses=nk, zRange=year.range, 
		winsize=window, n.rows=10, projSteps=100, 
		pr.rows=10, saveSize=10, numCores=ncores, saveDir=opath_it)

# use append.results to create single results file
	append.results(directory=opath_it, winsize=window, yearRange="2000-2017")

# use add.vars to create full variable set as individual result files
	data.file <- "ITvars_5win_2000-2017.rds"				# specify main results file to work from
	data.path <- paste0(opath_it, data.file)				# it should be in resultsDir
	add.vars(IT.file = data.path, ndvi.raster = meanndvi.path,
		class.raster = classes.path, filesDir = opath_it, 
		subDir.name = "ITvars")

# make rasters
	if (mk_ras == TRUE) {
		vars.path <- paste0(opath_it, "ITvars")	# specify folder holding files from add.vars

	# create dir for raster files.
		fold <- paste0(opath, "raster_results")
		subfold <- file.path(fold, "it")
		if(!dir.exists(fold)) {dir.create(fold)}
		if(!dir.exists(subfold)){dir.create(subfold)}

	IT.rasters(filesDir = vars.path, saveDir = subfold, minProp = 0.8, 
		Nyears = length(year.range), n.rows = 10, winsize = window, 
		class.raster = classes.path, raster.format="GTiff")
	}

