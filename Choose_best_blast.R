
# 2019-9-19 JCF

# R script to choose best hit from a blast tabular output (6) based on length of hit, % identity, 
# and bit score.

# Requires library GenomicsRanges
 library("GenomicRanges")

# Read in arguments
args <- commandArgs( trailingOnly=TRUE )
blast.file <- args[1]
prot.len <- args[2] * 3

# Read data in, force numeric columns as numeric
blast <- read.table(blast.file, header=FALSE)
blast[, 3:12] <- sapply(blast[, 3:12] , as.numeric)

# Check direction on reference
# Do matches go in same direction on ref?
blast <- cbind(blast, V13=sign( blast[,10] - blast[,9] ) )

# Set minimum useful identity of matches
blast <- blast[ blast[,3]>=50,]
# Reset row names after subset
rownames( blast ) <- seq( length <- nrow( blast) )

# Set ranges
i=1
blast.subs <- blast[ blast[,2]==  unique( blast[,2] )[i] ,]
ranges   <- IRanges(blast.subs[,7], blast.subs[,8] ) 

# Thout GRanges should be able to distinguish between scaffolds, but does not.
#ranges   <- GRanges( seqnames <- blast[,2],
#		     ranges   <- IRanges(blast[,7], blast[,8] ) 
#		  )


# Search ranges against itself to find all overlaps
range.overlaps <- cbind( subjectHits( GenomicRanges::findOverlaps(query=ranges, subject=ranges, minoverlap=20, type="any") ) , 
		queryHits( GenomicRanges::findOverlaps(query=ranges, subject=ranges, minoverlap=20, type="any") ) 
	       )


# If overlaps exist (excluding self-matches)
if ( sum( range.overlaps[,1]!=range.overlaps[,2]) > 0 ) {
	# Remove self-matches
	overlaps_to_check <- range.overlaps[ range.overlaps[,1]!=range.overlaps[,2] ,]
	
	# Compare percent identity of the two matches
	# Get the index of the lower scoring row within the subset
	k = 1
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
	
	}
    
    

# Get hits
hits = unique( blast[,2] )

# Initialize matrix to store results
qual_metrics = matrix(NA, ncol=5, nrow=length(hits))

# Loop through hits, subset the data for each, and calculate quality metrics.
for ( i in 1:length(hits) ) {
	blast.subs <- blast[ blast[,2]==hits[i],]
	# When matches hit the same ref in different directions, keep only the better match ( by % ID).
	if ( abs( sum(blast.subs[,13])) != nrow(blast.subs)  ) {
		# Is the 1 set better than the -1 set?
		if ( mean( blast.subs[ blast.subs[,13]==1 ,][,3] ) >= mean( blast.subs[ blast.subs[,13]==-1 ,][,3] ) ) {
			# If true: change subset to 
			blast.subs <- blast.subs[ blast.subs[,13] == 1 ,] } else {
			blast.subs <- blast.subs[ blast.subs[,13] == -1 ,] } 
		
	}
	
	# Perc identical matches
	qual_metrics[i,1] <- mean( blast.subs[,3] )
	
	# Length of match 
	qual_metrics[i,2] <- sum( abs( blast.subs[,10] - blast.subs[,9] + 1 ) )
	
	# Total number mismatches
	qual_metrics[i,3] <- sum( blast.subs[,5] )
	
	# Total number gaps
	qual_metrics[i,4] <- sum( blast.subs[,6] )
	
	# Total bit score
	qual_metrics[i,5] <- sum( blast.subs[,12] )
}

# When there are many hits on the same scaffold, the length & bit score 
# of poor matches can be better than the best hit. Filter poor hits out.
# When filtered down to 1 row, matrix becomes a vector, reformat matrix to get around this.
qual_metrics <- qual_metrics[ qual_metrics[,1] > 70 ,]
qual_metrics <- matrix(qual_metrics, ncol=5)

# Will decide based on % identity, length of match, & bit score
# If only 1 match survived identity filter, return that match.
# If all 3 metrics agree, return the match.
# Otherwise return "Fail"
if ( nrow(qual_metrics) == 1 ) { toString( hits[1] ) } else if
( which.max(qual_metrics[,1]) == which.max(qual_metrics[,2]) & which.max(qual_metrics[,2]) == 
		which.max(qual_metrics[,5]) ) { toString( hits[which.max(qual_metrics[,2])] )
} else { print("Fail") }




min( blast[ blast[, 2]== hits[1] ,9:10] )
max( blast[ blast[, 2]== hits[1] ,9:10] )


paste( hits[1], ":", min( blast[ blast[, 2]== hits[1] ,9:10] ), "-",  max( blast[ blast[, 2]== hits[1] ,9:10] ), sep="")
      
     

      
      
