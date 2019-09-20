
# 2019-9-19 JCF

# R script to choose best hit from a blast tabular output (6) based on length of hit, % identity, 
# and bit score.

# Read in arguments
args <- commandArgs( trailingOnly=TRUE )
file <- args[1]

# Read data in, force numeric columns as numeric
blast = read.table(file, header=FALSE)
blast[, 3:12] <- sapply(blast[, 3:12] , as.numeric)

# Get hits
hits = unique( blast[,2] )

# Initialize matrix to store results
qual_metrics = matrix(NA, ncol=5, nrow=length(hits))

# Loop through hits, subset the data for each, and calculate quality metrics.
for ( i in 1:length(hits) ) {
	blast.subs <- blast[ blast[,2]==hits[i],]
	
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
qual_metrics <- qual_metrics[ qual_metrics[,1] > 70 ,]


# Will decide based on % identity, length of match, & bit score
# If only 1 match survived identity filter, return that match.
# If all 3 metrics agree, return the match.
# Otherwise return "Fail"
if ( nrow(qual_metrics) == 1 ) { print qual_metrics } elif ( which.max(qual_metrics[,1]) == which.max(qual_metrics[,2]) & which.max(qual_metrics[,2]) == which.max(qual_metrics[,5]) ) {
	toString( hits[which.max(qual_metrics[,2])] )
} else { print("Fail") }








