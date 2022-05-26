# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.0        2015-08-20         Simon Hajjar        initial development
#     0.1        2018-12-17        Drew Thompson               moved to WDL

"""
sketchyRegion is a file used to mark sketchy regions when filtering a .VCF file
Created by Simon Hajjar on August 20, 2015
"""

import pysam

block = 1000 # The distance corresponding to a different neighbourhood of reads
false_negative = 0.75

def sketchyRegion(bam, record):
    """
    sketchyRegion(bam, record) returns true if the given VCF record
    is located in a sketchy area in the given .BAM file
    """

    if bam == ".":
        return "No_Bam"

    # Open the .bam file and fetch the appropriate reads
    bamfile = pysam.AlignmentFile(bam, "rb")

    chr1 = record.CHROM
    chr2 = record.INFO["CHR2"]
    pos1 = int(record.POS)
    pos2 = int(record.INFO["END"])

    # SKETCH TEST ONE: DID COVERAGE INCREASE DRAMATICALY IN THIS AREA?
    countMid = bamfile.count(chr1, pos1 - block,       pos1 + block      )
    countLef = bamfile.count(chr1, pos1 - (3 * block), pos1 - block      )
    countRig = bamfile.count(chr1, pos1 + block,       pos1 + (3 * block))
    if countMid > false_negative * (countRig + countLef): # Is coverage here drastically higher than expected?
        bamfile.close()
        return "BP1_Coverage"

    # SKETCH TEST TWO: DID COVERAGE INCREASE DRAMATICALLY IN THE SECOND BREAKPOINT?
    countMid = bamfile.count(chr2, pos2 - block,       pos2 + block)
    countLef = bamfile.count(chr2, pos2 - (3 * block), pos2 - block)
    countRig = bamfile.count(chr2, pos2 + block,       pos2 + (3 * block))
    print("%d %d %d" % (countLef, countMid, countRig))
    if countMid > false_negative * (countRig + countLef): # Is coverage here drastically higher than expected?
        bamfile.close()
        return "BP2_Coverage"

    bamfile.close()
    return "?"
