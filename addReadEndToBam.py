#!/usr/bin/env python

# 2019-9-13
# Report right-most position of alignment on reference from a bam
# Prints with YS:i: tag in bam 

# Source: Written by dariober in a SeqAnswers.com thread
# http://seqanswers.com/forums/showthread.php?t=51162

import pysam
import sys

insam= sys.argv[1]
samfile = pysam.AlignmentFile(insam, "rb")
outfile = pysam.AlignmentFile("-", "wb", template=samfile)

for aln in samfile:
    ys= aln.reference_end
    if not ys:
        ys= -1
    aln.setTag('YS', ys)
    outfile.write(aln)

samfile.close()
outfile.close()
sys.exit()
