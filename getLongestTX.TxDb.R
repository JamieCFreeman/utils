# 2020-5-6 JCF

# Function to 1 longest transcript of all genes from TxDB S4 object.
# Includes all geneIDs, so if only one transcript, returns that transcript.

# Returns matrix where column 1 is "gene_id", column 2 is "tx_name"

# Only tested on Aedes TxDB created from "Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gtf",
# be cautious 

# Requires library "GenomicFeatures"


getLongestTX.TxDb <- function( txdb ) {
  
  # Function expects a TxDb object, if one is not provided stop execution & provide error.
  if( class( txdb ) !="TxDb" ) { stop("Function expects a TxDb as input.") }
  
  # Get length of transcript of all genes
  TX.lengths <- GenomicFeatures::transcriptLengths(TxDb.name )
  #     tx_id   tx_name    gene_id nexon tx_len
  # 1     1 AAEL012102-RB AAEL012102     6   2159
  # 2     2 AAEL012102-RA AAEL012102     5   2147
  # 3     3 AAEL021681-RA AAEL021681     3    288
  # 4     4 AAEL012109-RB AAEL012109     5   2505
  # 5     5 AAEL026799-RA AAEL026799     2    814
  
  # Get 
  # duplicated() gives true for values that are duplicate of previous values, take IDs of genes that 
  #   have a duplicated value, when >2 transcripts, this includes the rest too, so need to take unique
  mult.TX.list <- unique( TX.lengths$gene_id[ duplicated( TX.lengths$gene_id ) ] )
  
  # Loop over genes that have multiple transcripts
  mult.TX.list.long <- cbind("GENEID"=mult.TX.list, "longest.TX"=NA)
  for (i in 1:length(mult.TX.list) ) {
    
     # Subset based on current id & take row w max transcript length
     # Get name of longest transcript
     longest.tmp <- which.max( TX.lengths[ TX.lengths$gene_id == mult.TX.list[i] ,]$tx_len )
     mult.TX.list.long[i,2] <- TX.lengths[ TX.lengths$gene_id == mult.TX.list[i] ,]$tx_name[longest.tmp]
    
  }
  

  # Final table combines the single transcript genes, the longest transcript of the multi-transcript genes
  longestTXperGene <- rbind(
                             # Single-isoform genes
                             TX.lengths[ !(TX.lengths$gene_id %in% mult.TX.list) ,],
                             # Longest transcripts
                             TX.lengths[ TX.lengths$tx_name %in% mult.TX.list.long[,2] ,]
                           )
  
  # Reorder to align with original input, include only columns of gene_id & txname
  longestTXperGene <- longestTXperGene[ order( longestTXperGene$tx_id ),] 
  longestTXperGene <- cbind("gene_id"=longestTXperGene$gene_id, "tx_name"=longestTXperGene$tx_name )
  
  
  # Function sanity checks
  # Output should have number of rows equal to number of unique gene_ids in input
  sanity <- nrow(longestTXperGene) == length( unique( TX.lengths$gene_id ) ) &&
            # Number of gene_ids should be == total # transcripts, minus transcripts beyond the first for each gene, 
            #     + the number of genes that have multiple transcripts (b/c we want to keep one for multi-transcript genes)
            length( unique( TX.lengths$gene_id ) ) ==   
            nrow(TX.lengths) - sum( table(TX.lengths$gene_id)[ ( table(TX.lengths$gene_id) > 1 ) ] ) +
            nrow(table(TX.lengths$gene_id)[ ( table(TX.lengths$gene_id) > 1 ) ])
  
  # If sanity check is not TRUE something went wrong. Error out, don't return results.
  # If sanity check is true, return results
  
  if ( sanity == FALSE ) {
  stop("Something went wrong") } else { return( longestTXperGene) }
  
  
}








#
# TODO: check timing on two approaches- may matter when going over all transcripts

# duplicated() gives true for values that are duplicate of previous values, take IDs of genes that 
#   have a duplicated value, when >2 transcripts, this includes the rest too, so need to take unique

#mult.TX.list <- unique( TX.lengths$gene_id[ duplicated(TX.lengths$gene_id) ] )

# # table() gives a count of how many times a particular CYP.TXSTART$GENEID occurs
# #   Take geneID that occur >1 time
# 
# mult.TX.list2 <- rownames( table(TX.lengths$gene_id)[ ( table(TX.lengths$gene_id) > 1 ) ] )

