
# 2019-9-19 JCF

# R script to choose best hit from a blast tabular output (6) based on length of hit, % identity, 
# and bit score.
# Tested for porting annotations between genome versions, so only looking for very high homology matches.
# Other use cases likely require different match criteria.
# Intended to compare tblastn/blastx output (protein to nt or nt to protein), so combines mulitple hits to same 
# subject for an overall score.
# Depends on additional script "Remove_multi_ref_hits_tblastn.R" to remove overlapping hits to query for
# multiple matches in subject (Finds hits that overlap by at least 20 nt, and removes lower identity match).
# If only one subject has >70% average identity for all hits, that subject is returned.
# 	Else if one subject has combined the highest  % identity, length of match, & bit score, it is returned. 
#	Else "fail" is returned. This is purposely strict to help examine the  match criteria.
# Blast seems to sort output by the lowest e-value entry overall is first, then all other hits for that subject.

# Example input file before "Remove_multi_ref_hits_tblastn.R":
# XP_005182218.1  35      100.000 366     0       0       1       366     1328307 1327210 0.0     758
# XP_005182218.1  35      97.887  142     3       0       364     505     1325051 1324626 2.24e-88        301
# XP_005182218.1  35      26.471  544     334     14      10      502     6856218 6857804 6.91e-39        156
# XP_005182218.1  35      27.711  415     269     11      21      411     7024765 7023542 9.29e-38        131
# XP_005182218.1  35      29.213  89      63      0       414     502     7023470 7023204 9.29e-38        50.1
# XP_005182218.1  35      24.812  532     354     14      7       502     7007803 7009368 3.20e-35        145
# XP_005182218.1  35      24.938  405     278     10      24      411     7019038 7017851 6.15e-31        117
# XP_005182218.1  35      29.167  72      51      0       414     485     7017783 7017568 6.15e-31        41.2
# XP_005182218.1  35      27.049  366     248     10      59      412     7042041 7040965 6.13e-22        104

# Example input file after "Remove_multi_ref_hits_tblastn.R":
# XP_005182218.1  35      100     366     0       0       1       366     1328307 1327210 0       758     -1
# XP_005182218.1  35      97.887  142     3       0       364     505     1325051 1324626 2.24e-88        301     -1

# Output is name of best subject.

# > sessionInfo()
# R version 3.5.0 (2018-04-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Matrix products: default
# BLAS/LAPACK: /programs/OpenBLAS/lib/libopenblas_sandybridgep-r0.2.18.so

# locale:
# [1] C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# loaded via a namespace (and not attached):
# [1] compiler_3.5.0



# Read in arguments 
# Call script with bash, ex:
#    Rscript /workdir/jcf236/utils/Choose_best_blast.R ./prot_all_max8_filt/${PROT_ID}_top_hit.blast1.filt
args <- commandArgs( trailingOnly=TRUE )
blast.file <- args[1]
prot.length < args[2]

# Read data in, force numeric columns as numeric
blast <- read.table(blast.file, header=FALSE)
blast[, 3:12] <- sapply(blast[, 3:12] , as.numeric)

# Check direction on reference
# Do matches go in same direction on ref?
blast <- cbind(blast, V13=sign( blast[,10] - blast[,9] ) )

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
	
	# % identical matches
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
# If identity & bit score agree, then check length. As long as the best alignment by other metrics
# covers >= 95% of length, return that. (Seems like a few bp differences in some alignments are causing fail).
# Otherwise return "Fail"
if ( nrow(qual_metrics) == 1 ) { toString( hits[1]) 
} else if ( which.max(qual_metrics[,1]) == which.max(qual_metrics[,2]) & which.max(qual_metrics[,2]) == 
		which.max(qual_metrics[,5]) ) { toString( hits[which.max(qual_metrics[,2])] ) 
} else if ( which.max(qual_metrics[,1]) == which.max(qual_metrics[,5]) &
	        qual_metrics[which.max(qual_metrics[,1]),2] >= ( 0.95 * prot.length) ) {
	     toString( hits[which.max(qual_metrics[,1])] ) 
} else { print("Fail"); print(qual_metrics) }



