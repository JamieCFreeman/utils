
# 2019-9-22 JCF

# Takes 2 arguments: a blast output file & a threshold. Checks if any matches above threshold
# exist & returns TRUE or FALSE

# Read in arguments 
# Call script with bash, ex:
#    Rscript /workdir/jcf236/utils/Choose_best_blast.R ./prot_all_max8_filt/${PROT_ID}_top_hit.blast1.filt
args <- commandArgs( trailingOnly=TRUE )
blast.file <- args[1]
threshold <- args[2]

# Define function
BLAST_bool_threshold_match <- function(blast.file, threshold) {
	
	# Read data in, force numeric columns as numeric
	blast <- read.table(blast.file, header=FALSE)
	blast[, 3:12] <- sapply(blast[, 3:12] , as.numeric)

	# Do any of the hits have a percent identity of >threshold?
	result <- sum( blast[,3] >= threshold ) > 0

	# Return boolean as numeric 0=FALSE, 1=TRUE
	cat( as.numeric( result ) )
}

# Call function on arguments provided
BLAST_bool_threshold_match(blast.file, threshold)




