cat('Loading custom functions...')

################################################################################
#   GENERAL PURPOSE CUSTOM FUNCTIONS                                           #
################################################################################

### Normalizing and centering function. Standardizes to range that is
###  centered on zero (mean is zero) with at least one extreme value 
###  that is equal to -32767 or +32767.
###  (short/signed integer range: -32767:32767)
normalize_center_sint <- function(m){
	m2 <- m[] - mean(m[], na.rm=TRUE)
	m3 <- as.integer(m2 / (max(abs(m2)) / 32767))

	return(m3)
}


### Function for determining if number is whole
is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


### Function that returns sum of all files within a directory (mB)
dsz <- function(x) {
  sum(file.size(list.files(x, full.names=TRUE))) / 1024^2
}


### Function that converts character variables to numeric and sums them
c2ns <- function(x) {
  eval(parse(text=paste0('as.numeric(',
			 paste(ls(envir=.GlobalEnv, pattern=x),
			       collapse=')+as.numeric('), ')')))
}


### Function that finds n intervals with equal numbers of observations from x
intv_starts <- function (x, n)
{
    breaks <- ggplot2::cut_number(x,        # Find n (left-closed) bins from x
				  dig.lab=nchar(max(x))+1, n=n, right=FALSE)

    freq <- table(breaks)                          # Tabulate freq across bins

    starts <- gsub('\\[', '',                      # Parse break start points
		   gsub(',[0-9].*', '', levels(breaks)))

    starts <- as.integer(ceiling(                  # Convert to integers
				 as.numeric(starts)))

    return(list(starts, freq))
}


### Function that returns indices in random order (used for distributing
###  data to cores in an even fashion)
RandSampIx <- function(m, num_chunks) {
  sample(1:num_chunks, nrow(m), replace=TRUE)
}


### Function that returns the gigabytes of available memory
checkRAM <- function(os_env) {
  if (os_env == 'Linux') {
    suppressMessages(library(doMC)) # Shared memory parallelism for Linux
    output <- as.numeric(system( # Detect available RAM
                         "awk '/MemAvailable/ {print $2}' /proc/meminfo",
                         intern=TRUE)) # Convert from kB to gB
  } else if (os_env == 'Windows') { # If using Windows start parallellism now
    suppressMessages(library(doSNOW)) # Shared memory parallelism for Windows
    x <- system2("wmic", args =  "OS get FreePhysicalMemory /Value",
                 stdout = TRUE) # Detect available RAM
    x <- x[grepl("FreePhysicalMemory", x)]
    x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
    x <- gsub("\r", "", x, fixed = TRUE)
    output <- as.numeric(x) # Convert from kB to gB
  }
}


### Check for package, install if missing, then load them
#list.of.packages <- c("ggplot2", "Rcpp")
checkLoad <- function(packs) {
  pack2inst <- packs[!(packs %in% installed.packages()[,"Package"])]
  if(length(pack2inst)) {
    print(paste('Installing', paste0(pack2inst, collapse=', ')))
    install.packages(pack2inst)
  }
  cat('Loading packages...')
  for (I in packs) {
    suppressMessages(library(I,                    # Load packages
			     character.only=TRUE))
  }
}


### Wrapper function for starting parallelism (Linux or Windows)
start_par <- function() {
  if (os_env == 'Linux') {
    registerDoMC(ncores)     # Register CPU cores
    cl <- NULL               # Unused dummy variable if in Linux
  } else if (os_env == 'Windows') { # If using Windows start parallellism now
    cl <- parallel::makeCluster(ncores)   # Make paral socket cluster
    registerDoSNOW(cl)       # Register CPU cores
  }
  return(cl)
}


### Wrapper function for terminating parallelism (only necessary for Windows)
stop_par <- function() {
  if (os_env == 'Windows') { # If using Windows stop parallellism
    stopCluster(cl)
  }
}


### Extract bounding box extent for a single state or single county.
###  crs.ll and crs.laea variables must be loaded (see vars.R)
exbb <- function(state, county) {
  cnty.ll <- ggplot2::map_data('county') # Get state, county boundaries
  if (nchar(cnty) > 0) { # If both state and county are requested
    sub.ll <- cnty.ll[cnty.ll$region == state &
                      cnty.ll$subregion == county, ]
  } else { # If only state is requested
    sub.ll <- cnty.ll[cnty.ll$region == state, ]
  }
  coordinates(sub.ll) <- ~long+lat       # Convert to spatial data frame
  proj4string(sub.ll) <- crs.latlon      # Specify coordinate reference system

  sub.laea <- spTransform(sub.ll, crs.laea) # Transform to Lambert AzmEqArea
  sub.laea.dfr <- data.frame(X=coordinates(sub.laea)[,1],
                             Y=coordinates(sub.laea)[,2],
                             group=sub.laea$group) # Copy to data frame

  # Extract boundaries for given state, county
  output <- extent(min(sub.laea.dfr$X), max(sub.laea.dfr$X),
               min(sub.laea.dfr$Y), max(sub.laea.dfr$Y))

  return(output)
}


################################################################################
#   FUNCTIONS FOR IMPUTING OVER MISSING VALUES                                 #
################################################################################

### Median value function for arrays
###  (Returns the middle number for sets with odd numbers of elements
###  (50th percentile), and for sets with even numbers of elements the value
###  representing the top of the lower range is returned.
###  e.g., for the set [2,3] 2 will be returned.
###  For the set [2,3,4] 3 will be returned.
###  For the set [2,3,4,5] 3 will be returned.)
median_val <- function(x, NAval=NULL) {

	if (!missing(NAval)) { # If optional argument declared...
		    x[which(x == NAval)] <- NA # Replace NAvals with NA
	}

	x_s <- apply(x, MARGIN=2, FUN=sort, na.last=FALSE) # Sort each column

	nvals <- function(x_s) {length(which(!is.na(x_s)))} # Fxn:Num nonNA vals

        y1 <- apply(x_s, MARGIN=2, FUN=nvals ) # Num nonNA vals in each column
        y2 <- y1 / 2 # Half y1
        mdr <- ceiling(nrow(x) - y2) # Rounded-up row no. of mdn val in ech col

	idx <- cbind(mdr, seq_along(mdr)) # 2 col matrix of relevant rows, cols

        output <- x_s[ idx ] # Middle value for each col
        output[which(y1 == 1)] <- NA # Set any col with only 1 nonNA val to NA

	return(output) # Return output
}


### Function for imputing where each missing value is replaced with its
### corresponding value from the median of all other years.
impute_median <- function(input, spc, Zth){
  # input: Input data array of one pixel per row
  # Zth: Threshold proportion of each pixel year that is allowed to
  #      contain zeroes without being imputed
  # spc: Samples per cycle (year)
  # OUTPUT VARIABLES (returned as a list)
  # output: Imputed (if possible) output with same dimensions as input
  # nImpVals: Number of values (elements) gap-filled in a previous
  #     iterations. (Fxn sums nImpVals from prior iteration with
  #     this iteration)

  nImpVals=0                                   # Default num imputed vals

  NAr <- mwhich(input, 1:ncol(input), NA,          # Pixels w/NAs
                'eq', 'OR')
  Zrsum <- rowSums(input == 0, na.rm=TRUE)         # Num of zeros for each pix
  Zr <- which((Zrsum >= spc)==T)                   # Pix with >= 1 yr of zeros
  NAZr <-  unique(c(NAr,Zr))                       # Pixels to go thru imptng 

  # If NA or excessive zeros exist then gap-fill
  if (any(length(NAZr))) {

    # Loop over pix requiring gap-filling
    for (I in 1:length(NAZr)) {

      # Reorder from 1 pixel/row to 1 pix-yr/row
      input2 <- matrix(t(input[NAZr[I], ]),        # Reformed matrix
                       ncol=spc, byrow=TRUE)

      Zrsumi3 <- rowSums(input2 == 0,na.rm=TRUE)   # Num zeros in each pix-yr

      Zri3 <- which((Zrsumi3 > spc*Zth) == TRUE)   # Rows > 3/4 zeros

      # If pix-yr is perdominantly zero-data then...
      if (length(Zri3) > 0) {
        input2[Zri3,] <- NA                        # ...Set those rows to NA
      }

      md <- median_val(input2, NAval=0)            # Calc col medians

      # If a complete year of median values can be made then gap-fill
      if (!anyNA(md)) {                            # If no NA's in md
        w  <- which(is.na(input2),                 # Index of NA's
                    arr.ind=TRUE)

        input2[w] <- md[w[,2]]                     # Replace NA element(s)
                                                   #  w/ corresponding median

        # Search for consecutive zeros & replace with respective median
        input2v <- as.vector(t(input2))            # Make temporary vector
        rlr <- rle(input2v)                        # Run len rpt cnt & values
        input2v <- as.vector(t(input2))            # Find consecutive repeats
        rpti <- which(rlr$lengths > spc/2          # Idx of zeros repeated
                      & rlr$values == 0)           #  more than 23 times
        if (length(rpti) > 0) {                    # If any consctv repeats..
          for (Z in 1:length(rpti)) {
            rptb <- sum(rlr$lengths[1:(rpti[Z]-1)])+1 # Beg idx of rpt vals
            rpte <- rptb + rlr$lengths[rpti[Z]] - 1   # End idx of rpt vals
            rptm <- rptb:rpte %% spc               # Modulus of indices
            rptm[rptm == 0] <- spc
            input2v[rptb:rpte] <- md[rptm]
          }
        }

        input[NAZr[I], ] <- input2v                # Update with imputed vals

        rm(input2v, rlr, rpti, rptbeg, rptm)

        nImpVals <- length(w)                  # Num of imputed vals
      }
      rm(list=Filter(exists, c('input2','md','w')))
    }                                              # End I loop
  }

  output <- list(imputed = input,                  # Arr w/ same dims as input
		 nImpVals = nImpVals)      # Num imputed values

  return(output)
}


impute_median_parallel <- function(input, ncores, spc, Zth){
  # input: Input data array of one pixel per row
  # Zth: Proportion of each pixel year that is allowed to contain zeros
  #      without being gap-filled
  # ncores: Number of processesors to use
  # spc: Samples per cycle (year)
  # OUTPUT VARIABLES (returned as a list)
  # output: Imputed (if possible) output with same dimensions as input
  # nImpVals: Number of values (elements) gap-filled in a previous
  #     iterations. (Fxn sums nImpVals from prior iteration with
  #     this iteration)
  # input_groups # Randomized group assignments for efficient parallelization

  # For efficiency randomize which cells are
  # assigned to which threads
  nImpVals=0                                   # Default num imputed vals
  input_groups <- RandSampIx (input, ncores)       # Randomized groups

  # Loop in parallel over input_groups
  unsorted <- foreach (P1 = 1:ncores, .packages=c('bigmemory'), .combine='rbind', .inorder=FALSE) %dopar% {
    # Make subset for this parallel thread
    input2 <- input[which(input_groups == P1), ]

    NAr <- mwhich(input2, 1:ncol(input2), NA,      # Pixels w/NAs
                  'eq', 'OR')
    Zrsum <- rowSums(input2 == 0, na.rm=TRUE)      # Num of zeros for each pix
    Zr <- which((Zrsum >= spc)==T)                 # Pix with >= 1 yr of zeros
    NAZr <-  unique(c(NAr,Zr))                     # Pixels to go through GF 

    # If NA or excessive zeros exist then gap-fill
    if (any(length(NAZr))) {

      # Loop over pix requiring gap-filling
      for (I in 1:length(NAZr)) {

        # Reorder from 1 pixel/row to 1 pix-yr/row
        input3 <- matrix(t(input2[NAZr[I], ]),     # Reformed matrix
                         ncol=spc, byrow=TRUE)

        Zrsumi3 <- rowSums(input3 == 0,na.rm=TRUE) # Num zeros in each pix-yr

        Zri3 <- which((Zrsumi3 > spc*Zth) == TRUE) # Rows > 3/4 zeros

        # If pix-yr is perdominantly zero-data then...
        if (length(Zri3) > 0) {
          input3[Zri3,] <- NA                      # ...Set those rows to NA
        }

        md <- median_val(input3, NAval=0)          # Calc col medians

        # If a complete year of median values can be made then gap-fill
        if (!anyNA(md)) {                          # If no NA's in md
          w  <- which(is.na(input3),               # Index of NA's
                      arr.ind=TRUE)

          input3[w] <- md[w[,2]]                   # Replace NA element(s)
                                                   #  w/ corresponding median

          # Search for consecutive zeros & replace with respective median
          input3v <- as.vector(t(input3))          # Make temporary vector
          rlr <- rle(input3v)                      # Run len rpt cnt & values
          input3v <- as.vector(t(input3))          # Find consecutive repeats
          rpti <- which(rlr$lengths > spc/2        # Idx of zeros repeated
                        & rlr$values == 0)         #  more than 23 times
          if (length(rpti) > 0) {                  # If any consctv repeats..
            for (Z in 1:length(rpti)) {
              rptb <- sum(rlr$lengths[1:(rpti[Z]-1)])+1 # Beg idx of rpt vals
              rpte <- rptb + rlr$lengths[rpti[Z]] - 1   # End idx of rpt vals
              rptm <- rptb:rpte %% spc             # Modulus of indices
              rptm[rptm == 0] <- spc
              input3v[rptb:rpte] <- md[rptm]
            }
          }

          input2[NAZr[I], ] <- input3v            # Update input2

          rm(input3v, rlr, rpti, rptbeg, rptm)

          nImpVals <- length(w)                # Num of vals gap-filled
        }
        rm(list=Filter(exists, c('input3','md','w')))
      }                                            # End I loop
    }

    input2 <- cbind(rep(P1, nrow(input2)),         # Compute last to guarantee
                    input2)                        #  shows up in foreach otpt
  }                                                # End foreach loop (P1)

  # Replace original input elements with gap-filled
  #  elements in their original order
  imputed <- unsorted[ , -1]                        # Remove a column
  imputed[] <- NA                                     # Initialize to NA
  for (J in 1:ncores) {
    imputed[which(input_groups == J), ] <-
            unsorted[which(unsorted[,1] == J), 2:ncol(unsorted)]
  }

  output <- list(imputed = imputed,
		 nImpVals = nImpVals)

  return(output)                                   # Return array of gap
                                                   # filled vals with same
                                                   # dims of initial input
}


################################################################################
#   FUNCTIONS FOR TRANSFORMING OUTPUT INTO RASTERS                             #
################################################################################

#####
#
# make.raster is a function for creating raster-format results:
# (does not save a new raster automatically).
# only works under assumption that 'data' has same length as !is.na in 'template', which it always should.
#####

# newdata: vector of cell values output raster will have; must have same length as non-NA cells in template.
# template: raster template, usually from raster.template function.
# varname: internally stored variable name.

make.raster <- function(newdata, template, varname) {
	# get template
		r.tmpl <- raster(template)	

        # get template vals and set to new vals
                rdat <- getValues(r.tmpl)
                rdat[!is.na(rdat)] <- newdata
                newr <- setValues(r.tmpl, rdat)
		names(newr) <- varname
	return(newr)
}


annual.rasters <- function(annual.data, data.name=NULL, yrs, template, 
		rtmpdir, remove.hrs, raster.format,
		save.dir=getwd()) {
	
	rasterOptions(tmpdir=rtmpdir)

	ptm0 <- proc.time()

	# create dir for these result files.
		fold <- file.path(save.dir, "raster_results")
		subfold <- file.path(save.dir, "raster_results", data.name)
		if(!dir.exists(fold)) {dir.create(fold)}
		if(!dir.exists(subfold)){dir.create(subfold)}

	dt <- readRDS(annual.data)
    	dt <- setDT(list(dt))
	rtempl <- raster(template)		# reference raster template
	
	if(names(dt)=="V1"){
		pfx <- sub(".rds$","",sub(".*/","",annual.data))	
		names(dt) <- pfx}
	dt[, yidx := rep(1:yrs, times=nrow(dt)/yrs)]
	pb = txtProgressBar(min = 0, max = length(yrs), initial = 0)
	for(i in 1:yrs){
		dat <- dt[yidx==i]
			for(j in 1:(ncol(dat)-1)){
				# make new raster object
					r <- make.raster(newdata=dat[[j]], template=template, varname=names(dat)[j])
				# write new raster with appropriate name and format.
					if(i<10){ rname <- paste0(names(dat)[j],".y0",i)}	# name by var and year
					if(i>9){ rname <- paste0(names(dat)[j],".y",i)}
					if(raster.format == "GTiff"){fext <- ".tif"}		# set extension for raster.format
					if(raster.format == "HFA"){fext <- ".img"}	
					if(raster.format == "raster"){fext <- ".grd"}	
					r.name <- paste0(rname, fext)
					r.path <- file.path(subfold, r.name)
					writeRaster(r, filename=r.path, format=raster.format, datatype='FLT4S', overwrite=T)
				}	
			setTxtProgressBar(pb,i)
			removeTmpFiles(h=remove.hrs)
		}
	rm(dt);rm(dat);rm(rtempl);rm(r);rm(rname);rm(r.name);rm(r.path)
	invisible(gc())
	ptm1 <- proc.time() - ptm0
}

###############

### annual.rasters.m saves annual raster output similar to annual.rasters, but where input
### data columns are read from a matrix (e.g., big data matrix), not from an rds file.
### one difference is the need here to specify column names (var.names).

annual.rasters.m <- function(annual.data, data.name=NULL, var.names, yrs,
			     template, datatype, rtmpdir,
			     remove.hrs, save.dir=getwd()) {
	
	rasterOptions(tmpdir=rtmpdir)

	nvars <- length(var.names)
	nyrs <- length(yrs)
	pb = txtProgressBar(min = 0, max = nvars * nyrs,
			    style = 3)

	# create dir for these result files.
		fold <- file.path(save.dir, "raster_results")
		subfold <- file.path(save.dir, "raster_results", data.name)
		if(!dir.exists(fold)) {dir.create(fold)}
		if(!dir.exists(subfold)){dir.create(subfold)}

	for(K in 1:nvars){

    	dt <- setDT(list(annual.data[,K]))
	names(dt) <- var.names[K]	
	rtempl <- raster(template)		# reference raster template
	dt[, yidx := rep(1:yrs, times=nrow(dt)/yrs)]

	for(i in 1:yrs){
		dat <- dt[yidx==i]
		# make new raster object
			r <- make.raster(newdata=dat[[1]], template=template,
					 varname=names(dat)[1])
		# write new raster with appropriate name and format.
			if(i<10){ rname <- paste0(names(dat)[1],"_y0",i)}	# name by var and year
			if(i>9){ rname <- paste0(names(dat)[1],"_y",i)}
			r.path <- file.path(subfold, rname)
			writeRaster(r, filename=r.path, format='GTiff',
				    datatype=datatype, overwrite=T)
			cnt <- (K - 1) * nyrs + i
			setTxtProgressBar(pb, cnt)
		}
	removeTmpFiles(h=remove.hrs)
	}
	close(pb)
	rm(dt, dat, rtempl, r, rname, r.path)
}
###############



###########	function cluster.mean, to get expected ndvi values for phenoclasses ##############
#		(or get expected value for any other metric)

	# pclasses: file path to cluster assignments (phenoclasses)
	# metric.vals: file path to metric (usually ndvi) values.

cluster.mean <- function(pclasses, metric.vals, vals.from.file=FALSE){

	### load data
		pc <- readRDS(pclasses)		# phenoclasses (clusters) by year
	# NDVI values for each year come from PolarMetrics.
		if(vals.from.file==TRUE){v <- readRDS(metric.vals)
		} else {v <- metric.vals}
		if(max(v) < 10001 & max(v) > 500){v <- v/100}	# rescale if values are large
	### get mean ndvi values by phenoclass, across all years.
		dt <- data.table(pc=pc,v=v)			
		dt.means <- dt[, .(vmean = mean(v)), by=pc][order(pc)]	# get means and sort
		
		return(dt.means)
}




################################################################################
#   INFORMATION THEORETIC FUNCTIONS                                            #
################################################################################
###################### Here are the large IT functions.

############ IT.priors function to calculate prior transition probabilities ###################

IT.priors <- function(rasterData, numClasses=500, zRange=NULL, numCores=4, pr.rows=500, saveDir=getwd()){
	ptm0 <- proc.time()
	cat("Calculating transition prior probabilities...",'\n')

	# set data source and make empty priors table
		rdat <- brick(rasterData)								# look at whole data set to count phenoclasses
		priors <- data.table(expand.grid(from=1:numClasses, to=1:numClasses))	# empty expansion table

	# set number of t-to-t+1 steps (year-to-year steps)
		if(is.null(zRange)) {yRange <- 1:(nlayers(rdat)-1)			# steps in foreach loop
			} else  yRange <- zRange[-length(zRange)]			

	# function to get transition frequencies
		transfreqs <- function(DATA){							# internal function to get freqs
			vals <- data.table(getValues(DATA))					# extract values and give columns convenient names
			names(vals) <- c("from","to")
			vals <- vals[complete.cases(vals),]					# remove cases with NAs
			frq <- vals[, .N ,by = list(from,to)]				# frequency table for unique transition types
			names(frq)[3] <- paste0("fq",i,".",k)
			rm(vals)
			return(frq)
		}

	# latitudinal sequences for raster row cropping 
		fromrow <- seq(1, nrow(rdat), pr.rows)			# raster row startpoints for data extraction
		torow <- fromrow + pr.rows - 1					# row endpoints
		torow[length(torow)] <- nrow(rdat)

	# do freq calculations for raster in latitudinal sequence
		pb = txtProgressBar(min = 0, max = length(fromrow), initial = 0) 
		cl=makeCluster(numCores)  								# set up for parallel
		registerDoParallel(cl) 									# register CPU cores
			for(i in 1:length(fromrow)){
				subr <- crop(rdat, extent(rdat, fromrow[i], torow[i], 1, ncol(rdat)))	# subset raster rows
					tab.freqs <- foreach(k = yRange, 						# loop
						.packages=c('data.table','raster')) %dopar%	{
							transfreqs(subr[[c(k,k+1)]])
							}
				for(j in 1:length(tab.freqs)){								# merge list results to 'priors'
					priors <- merge(priors, tab.freqs[j], by=c("from","to"), all.x=TRUE)
				}
				setTxtProgressBar(pb,i)
			}
		stopCluster(cl)											# end parallel

 	# sum freqs across all steps and drop step cols
		priors[is.na(priors)] <- 0													# replace NAs with zeroes
		priors[, freq := rowSums(.SD, na.rm = TRUE), .SDcols = 3:(ncol(priors))]	# sum relevant cols using data.table.
		priors <- priors[, .SD, .SDcols=c(1,2,ncol(priors))]						# drop step columns
		
	# probabilities from frequencies
		priors[, prob := freq/sum(freq)]											# as probabilities

	# save file
		filename <- paste0("transition_priors_periods",zRange[1],"-",zRange[length(zRange)],".rds")
		filename <- file.path(saveDir, filename)
		saveRDS(priors, filename)

	ptm1 <- proc.time() - ptm0
	cat('\n',"Finished calculating and saving priors.")
	cat('\n',round(as.numeric(ptm1[3]/60),2),"minutes to calculate priors",'\n')
	return(priors[,c(1,2,4)])
}

############ END IT.priors function ###########################################################



############ IT.do function to get IT values. #################################################

IT.do <- function(data, nsteps=projSteps, phen.vals=centroid.vals, priors=phen.priors){

		x <- data.table(data)

	# make basic from-to table from data
		cl1 <- unlist(as.list(x[,1:(ncol(x)-1)]))			# stack all to-from tranistions
		cl2 <- unlist(as.list(x[,2:ncol(x)]))
		fromto <- data.table(from=cl1, to=cl2)
		fromto <- fromto[!with(fromto,is.na(from)& is.na(to)),]
		y <- fromto[, .N ,by = list(from,to)]				# frequencies of unique transition types

	# get full matrix and do priors-filling.
		y[,pb := N/sum(N)]							# get probs from freqs for existing data
			f <- y[,unique(from)]						# find classes unique to one column or the other
			tu <- y[,unique(to)]
			f <- append(f,tu[tu%in%f==FALSE])				# append unique classes to the other column
			tu <- append(tu,f[f%in%tu==FALSE])	
		n <- data.table(expand.grid(from=f, to=tu))			# expand to all to-from combinations (all matrix cells)
		n <- merge(n, y, by=c("from","to"), all.x=T)			# get observed probs from y
		n <- merge(n, priors, by=c("from","to"), all.x=T)		# get priors
		n[is.na(pb), pb := prob]						# fill missing y vals with priors
		n[is.na(pb), pb := 1e-12]						# fill any remaining NAs with negligible pos val
		n[, pb := pb / sum(pb)]							# re-normalize probabilities

	# Make trans matrix (pb is probabilities, so this is the joint dstn)	
		jd <- xtabs(pb ~ from + to, data = n, na.action = na.omit)
	
	# matrix summaries
		rmd <- rowSums(jd, na.rm=T)						# Row marginal probs (froms)
		cmd <- t(as.matrix(colSums(jd, na.rm=T)))				# Column marginal probs (tos)
		c_r <- as.matrix(jd/rmd)   			              	# Conditional probs (col | row)

	# IT metrics
		ha <- sum(ifelse(cmd > 0, cmd*-log2(cmd), 0))			# shan entropy (H) for cmd
		hb <- sum(ifelse(rmd > 0, rmd*-log2(rmd), 0))			# shan entropy for rmd
		hx <- sqrt(ha*hb)								# geometric mean of ha and hb (meanH)
		mi <- sum(ifelse(jd > 0, jd*log2(jd / (rmd %*% cmd)), 0))	# mutual information (MI)
	
	# LDA and LDA-based meanH
		m.lda <- jd; diag(m.lda) <- 0		# LDA-only matrix, i.e., only off-diagonal values are non-zero.
		lda <- sum(m.lda)					# proportion of transitions that are offdiag. Landscape dynamic activity, LDA.
		if(lda > 0){
			jd.lda <- m.lda/lda	
			rmd.lda <- rowSums(jd.lda, na.rm=T)						# Row marginal probs (froms)
			cmd.lda <- colSums(jd.lda, na.rm=T)						# Column marginal probs (tos)
			ha.lda <- sum(ifelse(cmd.lda > 0, cmd.lda*-log2(cmd.lda), 0))	# shan entropy for cmd.lda
			hb.lda <- sum(ifelse(rmd.lda > 0, rmd.lda*-log2(rmd.lda), 0))	# shan entropy for rmd.lda
			h.lda <- sqrt(ha.lda*hb.lda)							# geometric mean of ha and hb: LDA meanH.
			}	else h.lda = 0
		h.lda.p <- h.lda * lda									# h.lda.p is 'Temporal H'.

	# Equilibrium projection and derived expected values of various metrics
		c_r.prj <- as.matrix(c_r %^% (nsteps+1))			# Project conditional probs n steps (years).
		rmd.prj <- as.numeric(rmd %*% c_r.prj)			# Project rmd by the prj matrix to get stable eq vector.
	
	# Compare rmd and projected equilibrium rmd (KL distance)
		kl.dist <- sum(ifelse(rmd.prj>0, ifelse(rmd>0, rmd.prj*log2(rmd.prj/rmd), 0),0))
	
	# H_ev and MI_ev	(expected eq values; can compare to meanH, MI)
		h.ev <- sum(ifelse(rmd.prj > 0, rmd.prj*-log2(rmd.prj), 0))
		mi.ev <- sum(ifelse(c_r > 0, (c_r*rmd.prj)*log2(c_r / rmd.prj), 0), na.rm=TRUE)	

	# mean current and expected equilibrium factor score and ndvi values.
		classes <- data.table(class=as.numeric(names(rmd)))										# all classes present in rmd/rmd.prj
		classes <- merge(classes,phen.vals, by.x="class", by.y="pc")					# factor score centroid vals for classes
		classes.rmd <- apply(classes[,2:ncol(classes)], 2, weighted.mean, w=as.vector(rmd))		# weighted mean factor score vals for rmd
		classes.prj <- apply(classes[,2:ncol(classes)], 2, weighted.mean, w=as.vector(rmd.prj))	# same for rmd.prj
	
	return(c(hx,mi,lda,h.lda,h.lda.p,kl.dist,h.ev,mi.ev,classes.rmd,classes.prj))				# return results.
}

########################################### END IT.do  ##################################



################### IT.parRun function to process data in parallel ######################

IT.parRun <- function(rasterData, numClasses=500, zRange=NULL, winsize=5, n.rows=5, projSteps=50,
		prior.probs=NULL, pr.rows=500, saveSize=100, numCores=detectCores()-1, saveDir=getwd()){
	
	# point to raster brick
		if(is.null(zRange)) {r <- stack(rasterData)}
		if(!is.null(zRange)) {r <- stack(rasterData, bands=zRange)}			# drop unwanted years

	# middle cell of window based on window size:
		midcell <- ceiling((winsize^2)/2)

	# sequences for raster row startpoint and file-save bookmarks
		starts <- seq(1, nrow(r), n.rows)				# raster row startpoints for data extraction
		nread <- c(rep(n.rows, length(starts)-1), nrow(r)-starts[length(starts)]+1)	# n.rows to read (final value varies)
		steps <- seq(1, length(starts), saveSize)		# start rows for result file save points

	# set up empty results table template.
		res.empty <- data.table(		
				hx=numeric(),mi=numeric(),lda=numeric(),hx.lda=numeric(),hx.lda.p=numeric(),
				kl.dist=numeric(),h.ev=numeric(),mi.ev=numeric(),f1=numeric(),
				f2=numeric(),f3=numeric(),f4=numeric(),ndvi=numeric(),f1ev=numeric(),
				f2ev=numeric(),f3ev=numeric(),f4ev=numeric(),ndvi.ev=numeric())

		# save empty template version for later use.
			saveRDS(res.empty, file=file.path(saveDir, "it.names.temp.rds"))

	# generate transition prior probabilities, using IT.priors (also saves an RDS file).
		if(is.null(prior.probs)){
			phen.priors <- IT.priors(rasterData=rasterData, numClasses=numClasses, 
			zRange=zRange, pr.rows=pr.rows, numCores=numCores, saveDir=saveDir)
		} else {phen.priors <- readRDS(prior.probs)}

	# loop large process steps, saving temp file each time
		ptm.all <- proc.time()
		cat("Calculating IT landscape metrics...",'\n')
		
		for (j in 1:length(steps)){
			ptm <- proc.time()						# time this step
			rowset <- starts[steps[j]:(steps[j]+saveSize-1)]	# define raster rows to get for step j
			rowset <- rowset[!is.na(rowset)]
			rowset.n <- nread[steps[j]:(steps[j]+saveSize-1)]
			rowset.n <- rowset.n[!is.na(rowset.n)]

			res.it <- res.empty						# start with empty results for step j

  		# Begin parallelism
			cl=makeCluster(numCores, type = "SOCK")
			clusterExport(cl, c("data.table","rbindlist","%^%",
			"projSteps","phen.priors","centroid.vals"), envir=environment())	

		# loop actual calculation process, executing IT.do function on sets of windows.
			pb = txtProgressBar(min = 0, max = length(rowset), initial = 0)
			for(i in 1:length(rowset)){
				vals <- getValuesFocal(r, row=rowset[i], nrows=rowset.n[i], 
						ngb=winsize, array=T)					# extract data in windows
				vals <- vals[which(!is.na(vals[,midcell,1])) ,, ]		# drop windows where focal cell is NA
						
				if(dim(vals)[1]>0){
					res.sub <- parApply(cl, vals, 1, IT.do)			# parallel apply function
					sub.vals <- data.table(t(res.sub))				# append results to results table
					l.vals <- list(res.it,sub.vals)
					res.it <- rbindlist(l.vals)
					}
				setTxtProgressBar(pb,i)
			}
			stopCluster(cl)				# end parallelism

			## save results 
			if(nrow(res.it) > 0){
				filename <- paste("ITvars_",winsize,"win_step",j,".rds",sep="")
				if(j < 10) {
					filename <- paste("ITvars_",winsize,"win_step0",j,".rds",sep="")
				}	
				filename <- file.path(saveDir, filename)
				saveRDS(res.it, file=filename)
				rm(res.it)
				gc()
			}
			ptm1 <- proc.time() - ptm
			cat('\n',as.numeric(ptm1[3]/60),"minutes for step",j,"of",length(steps),'\n')
		}
		cat('\n',"Finished calculating IT landscape metrics.",'\n')
		ptm.all1 <- proc.time() - ptm.all
		cat('\n',as.numeric(ptm.all1[3]/60/60),"hours to finish and save results.",'\n')
	## end processing loop. Nothing returned; files saved.
}
##################### END IT.parRun ################################################



################## append.results function to put all results into one file and delete temp result files ########

append.results <- function(directory, winsize, yearRange){
	cat("\n","Compiling results into single file...","\n")

	# get file list
		stepfiles <- list.files(path=directory, pattern = "_step", full.names=TRUE)
		
	# start empty results table, from template saved by IT.parRun.
		res.it <- readRDS(paste0(directory,"it.names.temp.rds"))
		
	# append files
		for (i in 1:length(stepfiles)){		# append the data
			sub.vals <- readRDS(stepfiles[i])
			l.vals <- list(res.it,sub.vals)
			res.it <- rbindlist(l.vals)
		}
	# save large file
		filename <- paste0(directory,"ITvars_",winsize,"win_",yearRange,".rds")
		saveRDS(res.it, file=filename)
	# delete step files
		if(file.access(filename)==0){unlink(stepfiles)}
	cat("\n","Finished compiling and saving single results file.","\n")
}
############################# END append.results ############################



###################### add.vars function to calculate derivative variables and save all single-var files #########

add.vars <- function(IT.file, ndvi.raster, class.raster, filesDir, subDir.name="ITvars"){
	
	ptm0 <- proc.time()
	cat('\n',"Calculating additional IT landscape metrics...",'\n')
	pb = txtProgressBar(min = 0, max = 40, initial = 0)

	dat <- readRDS(IT.file)							# result file from IT.parRun/append.results
	ndvi.r <- raster(ndvi.raster)					# scaling raster for ascendency etc.
	r.template <- raster(class.raster, band=1)		# raster mask

	# get matched ndvi values
		ndvi.r <- crop(ndvi.r, r.template)
		ndvi.m <- mask(ndvi.r, r.template)			# mask to cells with phenoclass data
		ndvi.vals <- getValues(ndvi.m)				# extract ndvi values as vector
		ndvi.vals <- ndvi.vals[!is.na(ndvi.vals)]	# remove NAs
			if(length(ndvi.vals)==nrow(dat))
				{dat[,meanndvi := ndvi.vals]} else	# add meanndvi to data
				{cat('\n',"Length of NDVI values doesn't match length of IT results",'\n')}

		# Keep in mind:
		# var 'ndvi' is estimated from current phenoclass composition;
		# var 'ndvi.ev' is estimated from equilibrium composition;
		# var 'meanndvi' is mean ndvi across years and pixels from the
		# original annual ndvi measures; it is used to scale IT vars.

	# add new vars and save them as individual rds files
	# (along the way, remove vars no longer needed)

		dir.create(paste0(filesDir,subDir.name))		# dir for all new files

		# define loop save function
			saveout <- function(vars){
				for(i in 1:length(vars)){	
					x <- dat[,.SD, .SDcols = (vars[i])]
					filename <- paste0(filesDir, subDir.name,"/", paste0(vars[i],'.rds'))
					saveRDS(x, file=filename)
				}
			}
	# save vars not used in further calcs
		savset <- c('lda','hx.lda','hx.lda.p','kl.dist')		# current working vars
			saveout(savset)							# use saveout func to save rds files	
			dat[, (savset) := NULL]						# done with these so remove
			setTxtProgressBar(pb,4)
	# factor score vars
		dat[,f1.eqdif := f1ev - f1]			# new vars based on factor scores
		dat[,f2.eqdif := f2ev - f2]
		dat[,f3.eqdif := f3ev - f3]
		dat[,f4.eqdif := f4ev - f4]
			savset <- c('f1','f2','f3','f4','f1ev','f2ev','f3ev','f4ev',
						'f1.eqdif','f2.eqdif','f3.eqdif','f4.eqdif')		# current working vars
			saveout(savset)
			dat[, (savset) := NULL]
			setTxtProgressBar(pb,16)
	# hx and me eqdifs
		dat[,hx.eqdif := h.ev - hx]
		dat[,mi.eqdif := mi.ev - mi]
			savset <- c('hx.eqdif','mi.eqdif')		# current working vars
			saveout(savset)
			dat[, (savset) := NULL]
			setTxtProgressBar(pb,18)
	# conditional entropy and overhead
		dat[,condent := hx - mi]
		dat[,condent.ev := h.ev - mi.ev]
		dat[,condent.eqdif := condent.ev - condent]
		dat[,ovhead := condent * meanndvi]
		dat[,ovhead.ev := condent.ev * ndvi.ev]
		dat[,ovhead.eqdif := ovhead.ev - (condent * ndvi)]
			savset <- c('condent','condent.ev','condent.eqdif','ovhead','ovhead.ev','ovhead.eqdif')
			saveout(savset)
			dat[, (savset) := NULL]
			setTxtProgressBar(pb,24)
	# development and capacity
		dat[,develop := ifelse(hx==0, 0, mi/hx)]
		dat[,develop.ev := ifelse(h.ev==0, 0, mi.ev / h.ev)]
		dat[,develop.eqdif := develop.ev - develop]
		dat[,capacity := hx * meanndvi]
		dat[,capacity.ev := h.ev * ndvi.ev]
		dat[,capacity.eqdif := capacity.ev - (hx * ndvi)]
			savset <- c('develop','develop.ev','develop.eqdif','capacity','capacity.ev','capacity.eqdif')
			saveout(savset)
			dat[, (savset) := NULL]
			setTxtProgressBar(pb,30)
	# ascendency and ndvi
		dat[,ascend := mi * meanndvi]
		dat[,ascend.ev := mi.ev * ndvi.ev]
		dat[,ascend.eqdif := ascend.ev - (mi * ndvi)]
		dat[,ndvi.eqdif := ndvi.ev - ndvi]
			savset <- c('ascend','ascend.ev','ascend.eqdif','ndvi','ndvi.ev','ndvi.eqdif')
			saveout(savset)
			dat[, (savset) := NULL]
			setTxtProgressBar(pb,36)
	# entropy and mutual information
		savset <- c('hx','mi','h.ev','mi.ev')
		saveout(savset)
		setTxtProgressBar(pb,40)
	ptm1 <- proc.time() - ptm0
	cat('\n',"finished")
	cat('\n',round(as.numeric(ptm1[3]/60),2),"minutes to calculate additional metrics and save files",'\n')
}
######################### END add.vars ########################################



################ IT.rasters function to make a raster file for each variable ###############
	
IT.rasters <- function(filesDir, saveDir, minProp=0, Nyears, n.rows=10, class.raster, 
				winsize, raster.format="GTiff", rtmpdir='F:/tmp'){
	
	ptm0 <- proc.time()
	cat('\n',"making and saving raster files for IT landscape metrics",'\n')
	### count transitions in each window, for masking.
	# function to get Ntransitions per window.
		IT.mask <- function(x){
			n.full <- sum(!is.na(x)) 			# count number of non-NA cells in window
			n.trans <- n.full*(Nyears-1)		# total transitions, not including NA cells
			return(n.trans)
		}
	
	# get raster template: class.raster
		r.template <- raster(class.raster, band=1)
	# specify center cell, given winsize.
		midcell <- ceiling((winsize^2)/2)
	# raster row startpoints for data extraction loop
		starts <- seq(1, nrow(r.template), n.rows)			
		nread <- c(rep(n.rows, length(starts)-1), nrow(r.template)-starts[length(starts)]+1)

	# empty results table.
		ntrans <- data.table(N=numeric())				

	# loop IT.mask function through all windows in raster.
		cat('\n',"Recording number of transitions in landscapes...",'\n')
		pb = txtProgressBar(min = 0, max = length(starts), initial = 0) 
		for(i in 1:length(starts)){
			# extract data in windows
				vals <- getValuesFocal(r.template, row=starts[i], nrows=nread[i], ngb=winsize)
				vals <- data.table(vals)
				vals <- vals[!is.na(vals[[midcell]])]					# remove NA windows
	
			if(dim(vals)[1]>0){
				res.sub <- apply(vals, 1, IT.mask)			# apply IT.mask
				sub.vals <- data.table(res.sub)				# append results to results table
				l.vals <- list(ntrans,sub.vals)
				ntrans <- rbindlist(l.vals)
			}
		setTxtProgressBar(pb,i)
		}
		cat('\n',"finished",'\n')
	### End counting transitions

	## calcs for masking by transitions counts (minProp specifies threshold)
		transmax <- ntrans[, max(N)]
		threshold <- minProp * transmax					# number of transitions needed for inclusion
		dat.pass <- ifelse(ntrans >= threshold, 1, 0)		# make vector indicating rows to keep
		prop.data <- sum(dat.pass)/nrow(ntrans)			# look at how much of the data remains
		cat('\n',round(prop.data*100,2), "% of pixels are within the defined threshold",'\n')

	### get data, create rasters, and save.
		cat('\n',"saving raster files",'\n')
		rasterOptions(tmpdir=rtmpdir)
	## variables data
		flist <- list.files(filesDir, pattern = "\\.rds$")
	## varnames for raster files
		flist.names <- sub('\\.rds.*', '', flist)			# strip file extension
	## feed data to raster template and save in specified file format.
		if(raster.format == "GTiff"){fext <- ".tif"}		# set extension for raster.format
		if(raster.format == "HFA"){fext <- ".img"}	
		if(raster.format == "raster"){fext <- ".grd"}	

		vals <- getValues(r.template)			# template raster values

		for(i in 1:length(flist)){
			dat <- readRDS(file.path(filesDir, flist[i]))					# get a variable
				if(is.data.frame(dat)==TRUE){dat <- dat[[1]]}
			dat.vals <- ifelse(dat.pass==1, dat, NA)						# mask var according to transmin
			vals.temp <- vals										# use copy of vals
			vals.temp[which(!is.na(vals.temp))] <- dat.vals					# set values for raster
			newrast <- setValues(r.template, vals.temp)					# feed values to raster
			rast.name <- paste(flist.names[i], "_",winsize,"win",fext,sep="")
			rast.name <- file.path(saveDir, rast.name)
			writeRaster(newrast, filename=rast.name, format=raster.format, overwrite=T)	#save raster
			removeTmpFiles(h=0.3)
		}
	ptm1 <- proc.time() - ptm0
	cat('\n',"finished")
	cat('\n',round(as.numeric(ptm1[3]/60),2),"minutes to make and save raster files",'\n')
}
############################ END IT.rasters #######################################################




########### function landscape.mean, to get metric mean across n by n pixel landscape in moving window #####
landscape.mean <- function(metric.vals, winsize, zRange=NULL, n.rows=10){
	cat('\n',"Loading data...",'\n')
		midcell <- ceiling(winsize^2/2)
		cellN <- winsize^2	
		cellNy <- cellN*length(zRange)

	### load data
		v <- stack(metric.vals[zRange])			# use zRange to stack layers you want.

	### extract appropriate values for windows and get means, with parallel processing.
		starts <- seq(1, nrow(v), n.rows)			# raster row startpoints for data extraction
		nread <- c(rep(n.rows, length(starts)-1), 	# n.rows to read (final value varies)
			nrow(v)-starts[length(starts)]+1)	

		# set up empty results table
		win.means <- data.table(winmean=numeric())
	cat('\n',"Calculating landscape mean values for IT scaling factor...",'\n')
	pb = txtProgressBar(min = 0, max = length(starts), initial = 0) 
			for(i in 1:length(starts)){
				#ptm0 <- proc.time()						# time this set			
				# extract data in windows
				vals <- getValuesFocal(v, row=starts[i], nrows=nread[i], ngb=winsize)
					vals <- setDT(vals)
					win.x <- vals[, .(winmean=rowSums(.SD, na.rm=T), 
						num.obs = ncol(vals) - Reduce("+", lapply(.SD, is.na)))]
					icell <- rep(1:(nrow(win.x) / cellN), cellN)
					win.x <- win.x[, lapply(.SD, sum), by=icell]
					win.x <- win.x[, .(winmean = winmean/num.obs)]
					l.vals <- list(win.means,win.x)				# append results to win.means
					win.means <- rbindlist(l.vals)
					setTxtProgressBar(pb,i)
				}
		win.means <- win.means/100				# rescale if data are in large-value format
	### write result as raster
		r.tmpl <- raster(v,1)							# get template
		newr <- setValues(r.tmpl, win.means[[1]])		# set to new vals
		newr <- mask(newr, r.tmpl)						# ensure NAs
		return(newr)
}
############# end landscape.mean ##########################################################################

cat('Done.\n')
