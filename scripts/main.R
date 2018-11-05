rm(list=ls())
invisible(gc())
print(date())

################################################################################

### Load program variables and custom functions
source('fxns.R')
source('vars.R')

### Print status messages
print(paste(rep('#', 73), collapse=''))
print(paste(rep('#', 73), collapse=''))
print(paste('Continue from prior run?:', continue))
print(paste('Restart codes:', paste(restart, collapse=' ')))
print(paste('gB of RAM available:',
	    sprintf('%0.1f', ram_avail_kb/1024^2)))
print(paste(date(), '| Clusters:', nk))
print(paste('Num seed trials:', iter_seed,
	'| Max num trials:', iter_trial,
	'| Num cores to use:',ncores))
print(paste(rep('#', 73), collapse=''))
print(paste('Mask water, ice, barren pixels using raster mask:', mask))
print(paste('Save output as rasters:', mk_ras))
print(paste(rep('#', 73), collapse=''))
print(paste('NDVI raster input directory:', srcpath))
print(paste('Output directory:', opath))
print('(Check here during preprocessing for % completion)')
print(paste(rep('#', 73), collapse=''))

t_tot0 <- proc.time()[3]                           # Mark time

################################################################################
### STEPS 1 - 3, PREPROCESS DATA: LOAD DATA, GAP-FILL, POLAR CALCULATION

print(paste(rep('#', 73), collapse=''))

source('preprocess.R')                             # Run preprocessing module
source('diag.R')                                   # Update diagnostics page

print(paste('Hrs gap-filling NDVI: ', helap_impute))
print(paste('Hrs calculating polar metrics: ', helap_calc_pm))

################################################################################
### STEP 4, CALCULATE FACTOR SCORES

print(paste(rep('#', 73), collapse=''))

source('fscore.R')                                 # Run FA, calc factor scores
source('diag.R')                                   # Update diagnostics page

print(paste('Done. Hrs calculating fact. scores: ', helap_calc_fs))

################################################################################
### STEP 5, CLUSTER FACTOR SCORES

print(paste(rep('#', 73), collapse=''))

source('cluster.R')                                # Run clustering module
source('diag.R')                                   # Update diagnostics page

print(paste('Done. Hrs clustering: ', helap_calc_fs))

################################################################################
### STEP 6, CALCULATE INFORMATION THEORETIC MEASURES

print(paste(rep('#', 73), collapse=''))

source('IT.R')
source('diag.R')                                   # Update diagnostics page

helap_tot <- (proc.time()[3] - t_tot0) / 60 ^2
print(paste('Tot runtime:', helap_tot))
