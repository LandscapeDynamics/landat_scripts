##############################################################################
### Gather statistics and update diagnostics html pages 
cat('Updating diagnostics...')

ncell <- readRDS(paste0(opath,             # Tot no. of pixels processed
			  'ncell.rds'))
nv <- as.double(ncell) * nyr * spc                # Tot no. of values processed

nIV <- sum(big.nImpVals[], na.rm=T)          # Tot. num of pix-yrs imputed

ts_data <- rep(0, nsamp-spc)

nv_ex <- (ncell - length(cid_good)) * nyr * spc    # No. of elements excluded

file <- readLines(paste0(dgpath,            # Read input
			       'bak/index.html'), n=-1)

file <- sub('Millions of Values:',   # Update text with num vals processed
	paste('Millions of Values:',
      sprintf('%0.1f', nv/10^6)), file)

if (restart$fs_complete == TRUE & nfac == 4) {
  load <- fa$loadings
  for (I in 1:6) {
	  tmp_out <- paste0(sprintf('%0.2f', f$loadings[I, ]), collapse='</TD><TD>')
	  tmp_out <- paste0('<TR><TD>',names(load[I]),'</TD><TD>', tmp_out, '</TD></TR>')
	  file <- sub(paste0('^<TR><TD>fv',I,'</TD></TR>'),
		      tmp_out, file)
  }
  writeLines(file, paste0(dgpath,           # Write output file
			  'index.html'))
}

cptrn <- c('pie_1_data$',                   # Character pattern to match
	   'bar_1_y$', 'bar_2_y$',          #  in HTML file
	   'bar_3_y$', 'bar_4_y$',
	   'ts_1_y$')

# Specify replacement text
#new <- 'data: [12, 19, 3, 5, 2, 3], // bar_1_y'
new_line <- vector(length=length(cptrn),     # Initialize variable
		   mode='character')

new_line[1] <- paste0('data: [',             # Insert character string
		      sprintf('%0.1f', nv/10^6), ', ',
		      sprintf('%0.1f', nIV/10^6), ', ',
		      sprintf('%0.1f', nv_ex/10^6), '], // pie_1_data')

new_line[2] <- paste0('data: [',             # Insert character string
		      helap_impute, ', ', helap_calc_pm, ', ',
		      helap_calc_fs, ', ', helap_calc_cls, ', ',
		      helap_calc_it, '], // bar_1_y')

new_line[3] <- paste0('data: [',             # Insert character string
#		      c2ns('^helap_load_'),', ',c2ns('^helap_wr_rds_'), ', ',
		      c2ns('^helap_wr_ras_'),', ',c2ns('^helap_pop_'), ', ',
		      c2ns('^helap_calc_'), '], // bar_2_y')

new_line[4] <- paste0('data: [',             # Insert character string
		      dsz(srcpath), ', ', dsz(opath_pm), ', ',
		      dsz(opath_fs), ', ', dsz(opath_cls), ', ',
		      dsz(opath_it), '], // bar_3_y')

new_line[5] <- paste0('data: [',             # Insert character string
		      0, ', ',
		      dsz(paste0(opath_ras,'pmetric/')), ', ',
		      dsz(paste0(opath_ras, 'fscore/')), ', ',
		      dsz(paste0(opath_ras, 'cluster/')), ', ',
		      dsz(paste0(opath_ras, 'it/')), '], // bar_4_y')

new_line[6] <- paste0('data: [',             # Insert character string
		      paste(ts_data, collapse=', '),
		      '], // ts_1_y')

for (I in 1:length(cptrn)) {
  idx <- grep(cptrn[I], file)               # Line no. to edit

  file[idx] <- new_line[I]                   # Update line

  writeLines(file, paste0(dgpath,           # Write output file
			  'index.html'))
}

### update run ID and date
fname <- c('index')
for (FN in fname) {
  tx <- readLines(paste0(dgpath, FN,'.html')) # Original HTML file
  tx2 <- sub(pattern = 'RUN_ID',                    # Replace text
             replace = run_id, x = tx)
  tx3 <- sub(pattern = 'DATE',                      # Replace text
             replace = date(), x = tx2)
  writeLines(tx3, con=paste0(dgpath, FN, '.html'))  # Write output file
}
rm(FN,tx,tx2,tx3)

#sum(file.size(list.files('../output/fs/', full.names=TRUE))) / 1024^3
cat('Done.\n')
