
# 2019-9-19 JCF

# From a set of space delimited ranges (generally the parsed output of a blast command),
# returns the length of the regions. For approx protein length (nt/3) set protein to TRUE

# Expected space-delimited input with interval partners within a row:
# 	24063005 24064420
#	24064479 24065201


matchLength <- function(file, protein=FALSE) {

	# Read in space-separated file.
	spans <- read.table(file, sep=" ")

	# Create vector to store 
	stor_vec <- rep(NA, nrow(spans) )

	# Get span length by row 
	for ( i in 1:nrow(spans) ) {
		stor_vec[i] <- abs( spans[i,1] - spans[i,2] + 1)
	}

	# If protein, want to divide by 3, if not don't
	if ( protein == TRUE ) {
		sum(stor_vec)/3
	} else { sum(stor_vec) }

}


