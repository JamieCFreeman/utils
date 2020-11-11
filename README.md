# utils

Useful functions for manipulating sequence data, mostly for my own reuse.
These were all written for specific use cases, so use at your own risk.



| Command | Description |
| --- | --- |
| BLAST_bool_threshold_match.R | Takes BLAST output (6) & checks if any matches >= % identity threshold, return boolean |
| Choose_best_blast.R  | Takes tblastn/blast output (6) & chooses the best one, gets overall score (across fragmented matches, eg exons)                            based on total % ID, bit score, & length of match  |
| Remove_multi_ref_hits_tblastn.R | From BLAST+ output for tblastn, find multiple hits of the same reference to the same sequence. For any regions with two overlapping hits, drop the lower scoring. |
| SRP2SRR_convert.R | query SRA to get a text file list of the SRR runs within that project |
| addReadEndToBam.py | Ammend BAM to include POS of rightmost alignment |
| check_inst.R | Check package install in R, return boolean |
| getLongestTX.TxDb.R | For transcript database (TxDB) S4 object in R, filter to get max 1 Tx per gene (picks longest) |
| my_colors.R | Source in other R functions for my personal color palettes  |

