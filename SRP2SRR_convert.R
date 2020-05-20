F
# 2019-6-20
#
# Using the R package SRAdb, query SRA to get a text file list of the SRR runs within that project
# Slow (must download whole SRA SQL table), but only 1X per project, so not a problem. 
# If need to use more frequently, maybe move to edirect utils (sometimes buggy for me, but faster)

# Dependencies:
	# R
	# SRAdb package, installed through Bioconductor if necessary
				# packageVersion('BiocManager')
				# BiocManager::install('SRAdb')
				
				
# Input:
	# BioProject accession number (SRP......)
	
# Output:
	# Txt file with one SRR entry per line.


# Using the SRAdb package in R

# Read in arguments 
# Call script with bash, ex:
# Rscript /workdir/jcf236/utils/SRP2SRR.R SRP015948
args <- commandArgs( trailingOnly=TRUE )
study.srp <- args[1]

# Load dependencies
library('SRAdb')

SRP2SRR <- function( study.srp ) {

	# Download the lite sql table from SRA to query
	# Get a list of the files contained in the project "SRP015948"
	sqlfile <- file.path(system.file('extdata', package='SRAdb'),'SRAmetadb_demo.sqlite')
	sqlfile <- getSRAdbFile()
	sra_con <- dbConnect( SQLite(), sqlfile )

	conversion <- sraConvert( study.srp, sra_con = sra_con )
	
	# Write run IDs to file
	write.table( x = conversion$run, file = paste(study.srp, ".txt" , sep = ""), row.names = FALSE, col.names = FALSE, 
					quote = FALSE)

}

SRP2SRR( study.srp )

# Deletion of SQL file. 
unlink( "SRAmetadb.sqlite" )


