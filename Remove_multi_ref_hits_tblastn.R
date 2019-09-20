
# 2019-9-20 JCF

# From BLAST+ output for tblastn, find multiple hits of the same reference to the same sequence.
# For any regions with two overlapping hits, drop the lower scoring.

# Requires library GenomicsRanges
library("GenomicRanges")

# Read in arguments
# args <- commandArgs( trailingOnly=TRUE )
# blast.file <- args[1]


blast.file="XP_011294256_top_hit.blast1"
# Read data in, force numeric columns as numeric
blast <- read.table(blast.file, header=FALSE)
blast[, 3:12] <- sapply(blast[, 3:12] , as.numeric)

# Do matches go in same direction on ref?
blast <- cbind(blast, V13=sign( blast[,10] - blast[,9] ) )

# Set minimum useful identity of matches
# blast <- blast[ blast[,3] >= 50,]
# Reset row names after subset
rownames( blast ) <- seq( length <- nrow( blast) )

# Set ranges
i=1
blast.subs <- blast[ blast[,2]==  unique( blast[,2] )[i] ,]
ranges   <- IRanges(blast.subs[,7], blast.subs[,8] ) 

# Search ranges against itself to find all overlaps
range.overlaps <- cbind( subjectHits( GenomicRanges::findOverlaps(query=ranges, subject=ranges, minoverlap=20, type="any") ) , 
		queryHits( GenomicRanges::findOverlaps(query=ranges, subject=ranges, minoverlap=20, type="any") ) 
	       )

# If overlaps exist (excluding self-matches)
if ( sum( range.overlaps[,1]!=range.overlaps[,2]) > 0 ) {
	# Remove self-matches
	overlaps_to_check <- range.overlaps[ range.overlaps[,1]!=range.overlaps[,2] ,]
	
	while ( nrow(overlaps_to_check) > 0 ) {
	
	# Compare percent identity of the two matches
	# Get the index of the lower scoring row within the subset
	k = 1
	# Each true match will be present in each direction.
	# Find the row that is the mirror of the current row and remove it.
	y <- rev( overlaps_to_check[k,] ) 
	mirror <- which( apply( overlaps_to_check, 1, function(x) identical( x , y) ) )
	overlaps_to_check <- matrix( overlaps_to_check[-mirror, ], ncol=2)
	
	subs.to.remove <- overlaps_to_check[k, which.min( 
						   c( blast.subs[ overlaps_to_check[k,1] ,3], 
						      blast.subs[ overlaps_to_check[k,2] ,3] 
						     ) 
						  )
				       ]

	# Get the index of the lower scoring row within the main set
	row.to.remove <- as.numeric( as.character( row.names( blast.subs[subs.to.remove, ] ) 
						 ) 
				   )

	# Remove the lower scoring row.
	blast <- blast[-row.to.remove, ]
	
	# Remove the fixed overlap from our list to check		       
	overlaps_to_check <- overlaps_to_check[-k,]

	}

}
