
# 2019-9-20 JCF

# From BLAST+ output for tblastn, find multiple hits of the same reference to the same sequence.
# For any regions with two overlapping hits, drop the lower scoring.

# rm( list = ls())
# Requires library GenomicsRanges
library("GenomicRanges")

# Read in arguments
args <- commandArgs( trailingOnly=TRUE )
blast.file <- args[1]
out.file <- args[2]

# blast.file="XP_005190501_top_hit.blast1"
# out.file="./recip_blast_out/XP_005190501_top_hit.blast1.filt"	
# out.file="./prot_all_max8_filt/XP_005182218_top_hit.blast1.filt"	# for blast1
# blast.file="XP_005182218_top_hit_pipeless.blast2"			#for blast2
#out.file="./recip_blast_out/XP_011294256_top_hit_pipeless.blast2.filt"	#for blast2
# Read data in, force numeric columns as numeric
blast <- read.table(blast.file, header=FALSE)
blast[, 3:12] <- sapply(blast[, 3:12] , as.numeric)

# Are there any good matches present? If not, break script.
# blast[,3] >= 50

# IRanges won't take negative ranges, add row to indicate sense (+1 for positive, -1 for negative)
blast <- cbind(blast, V13=sign( blast[,8] - blast[,7] ) )

# Set minimum useful identity of matches
blast <- blast[ blast[,3] >= 50,]
# Reset row names after subset
rownames( blast ) <- seq( length <- nrow( blast) )

# Set ranges 
for ( i in 1:length( unique( blast[,2] ) ) ) {
	blast.subs <- blast[ blast[,2]==  unique( blast[,2] )[i] ,]
	# Ensure ranges go in + direction
	if ( blast.subs[1,13] == 1) {
		ranges   <- IRanges(blast.subs[,7], blast.subs[,8] ) } else if ( blast.subs[1,13] == -1) {
		ranges   <- IRanges(blast.subs[,8], blast.subs[,7] )
		}

	# Search ranges against itself to find all overlaps
	range.overlaps <- cbind( subjectHits( GenomicRanges::findOverlaps(query=ranges, subject=ranges, minoverlap=20, type="any") ) , 
			queryHits( GenomicRanges::findOverlaps(query=ranges, subject=ranges, minoverlap=20, type="any") ) 
	       		)

	# If overlaps exist (excluding self-matches)
	if ( sum( range.overlaps[,1]!=range.overlaps[,2]) > 0 ) {
		# Remove self-matches
		overlaps_to_check <- range.overlaps[ range.overlaps[,1]!=range.overlaps[,2] ,]
		
		# While there are still rows in overlaps_to_check
		while ( nrow(overlaps_to_check) > 0 ) {
	
			# Initialize k 
			k = 1
			# Each true match will be present in each direction.
			# Find the row that is the mirror of the current row and remove it.
			y <- rev( overlaps_to_check[k,] ) 
			mirror <- which( apply( overlaps_to_check, 1, function(x) identical( x , y) ) )
			overlaps_to_check <- matrix( overlaps_to_check[-mirror, ], ncol=2)
					       
			# Compare percent identity of the two matches
			# Get the index of the lower scoring row within the subset
			subs.to.remove <- overlaps_to_check[k, which.min( 
						   c( blast.subs[ overlaps_to_check[k,1] ,3], 
						      blast.subs[ overlaps_to_check[k,2] ,3] 
						 	    ) 
							  )
				    	   ]
					   
			# Get the index of the lower scoring row within the main set
			row.to.remove <- row.names( blast.subs[subs.to.remove, ] ) 
							

			# Remove the lower scoring row, and resubset blast
			blast <- blast[ row.names(blast) != row.to.remove ,]

			subs.blast <- blast[ blast[,2]==  unique( blast[,2] )[i] ,]
	
			# Remove the fixed overlap from our list to check ( remove all rows containing it)
			# Must be last action- if nrow==0, loop breaks
			overlaps_to_check <- overlaps_to_check[ which(  !apply( overlaps_to_check, 1, 
						   function(x) identical( x[1] , subs.to.remove ) | identical( x[2] , subs.to.remove ) )
							) ,]
			# overlaps_to_check <- overlaps_to_check[-k,]

		}
	}

}
					       
# Write file to output
write.table(blast, file=out.file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)	       
					       
					       
					       
