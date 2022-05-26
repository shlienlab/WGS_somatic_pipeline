#!/usr/bin/python
# Written using Python v3.4.0

# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.0        2015-05-25         Simon Hajjar        initial development
#     0.1        2018-12-17        Drew Thompson               moved to WDL
#     0.2        2019-07-31        Drew Thompson      remove harcoded paths
#     0.2        2019-08-22        Drew Thompson Save filter info for graph

# vcf2Tab.py is an executable python file used to convert a .vcf file to a .tab
# file, adding the name of the genes associated with each column. Includes any
# additional fields added to the .vcf file, which may also contain multiple samples.

# Written by Simon Hajjar on 25 May, 2015

import sys
import os
import vcf
import csv
import numpy as np
from optparse import OptionParser
from sketchy import sketchyRegion
from UCSC_naming import geneName, geneBetween, exonName, exonBetween, polymorphisms

NULL_VAL = "0"    # Default used when no value is found for a cell
LIST_CONCAT = ";" # The string used to concatenate lists in a column
MAX_LIST = 10     # The maximum amount of genes allowed in-between 

COLUMNS = [["Chromosome1", "Chromosome2", "Position1", "Position2", "Breakpoint1", 
            "Breakpoint2",  "BAF_1", "Break_1", "Comments_1", "BAF_2", "Break_2", 
            "Comments_2", "Between", "Exons", "Common", "PE", "SR", "SVLEN", "ID", 
            "Reference", "Alternate", "Quality", "Filter", "Sketchy?", "INSLEN", "MAPQ", 
            "SVTYPE", "SVMETHOD", "PRECISE", "IMPRECISE", "CIPOS", "CIEND", "CT", "SRQ", "SAMPLE1_RC",
            "SAMPLE1_RR", "SAMPLE1_DV", "SAMPLE1_GQ", "SAMPLE1_RCR", "SAMPLE1_CN", "SAMPLE1_GT",
            "SAMPLE1_RCL", "SAMPLE1_DR", "SAMPLE1_FT", "SAMPLE1_GL", "SAMPLE1_RV", 
            "CONSENSUS"]]

def vcf2Tab(file, lowqual, bam, output_file, ucsc_info_path, filter_info_file, sv_type):
    """
    vcf2Tab(file, output, name) takes in a .vcf file and processes it
    into a .tab file, with gene names added to the columns.

    Args:
        vcf (string): the location of the .vcf file
        output (string): the directory to store the output in
        name (string): the name of the file to store results in
        bam (string) : the path to the .bam file
        lowqual (bool) : filter stuff out based on low quality

    Effects:
    creates a .tab file that is processed from the .vcf
    """

    print("Converting .vcf file to .tab file...")
    allowed = [str(i) for i in range(1, 23)] + ["X", "Y"]
    allowed += ["chr" + str(i) for i in range(1, 23)] + ["X", "Y"]

    # Write up the first 9 standard columns for the .tab file    
    tabFile = COLUMNS
    vcf_reader = vcf.Reader(open(file, 'r'))
    unfiltered_count = 0
    qual_filtered_count = 0
    for record in vcf_reader:

        if not record.CHROM in allowed:
            continue

        unfiltered_count += 1

        # Filter low-quality reads
        if lowqual:
            if (((record.INFO["PE"] < 2) and (record.INFO["SR"] < 2)) 
               or (record.INFO["MAPQ"] < 20)): continue
        else: 
            if "LowQual" in record.FILTER: continue

        qual_filtered_count += 1

        # Is it a sketchy region?
        sketchy = "?"
        try:
            sketchy = sketchyRegion(bam, record)
        except:
            sketchy = "Error"

        row = []

        # Adds any dictionary to the row, making sure to align correctly.
        # Helper function declared so that as many samples as the user likes
        # may be used in the final .tab file.
        def addDictToTab(dict):
            for field in dict:
                # If the field is not already in the table, we add a column
                pos = NULL_VAL
                if field == 'SR\n':
                    field = 'SR'
                elif field == 'CHR2':
                    pos = tabFile[0].index("Chromosome2")
                elif field == 'END':
                    pos = tabFile[0].index("Position2")
                #elif not field in tabFile[0]:
                #    tabFile[0].append(field)
                if pos == NULL_VAL:
                    pos = tabFile[0].index(field)

                # Extends the row in question to match up with the new add
                while len(row) <= pos:
                    row.append(NULL_VAL)

                # Inserts the given dictionary entry for this row. If the
                # entry is None, we use the NULL_VAL. If it is a list, we
                # concatenate them. Else, just put the value in the celll
                to_Insert = dict[field]
                if to_Insert == []:
                    to_Insert = NULL_VAL
                if type(to_Insert) is list:
                    to_Insert = LIST_CONCAT.join([str(i) for i in to_Insert])
                if to_Insert is None:
                    to_Insert = NULL_VAL
                row[pos] = to_Insert

        # Add the record information from the .vcf into the table
        addDictToTab(record.INFO)

        # Add the other information to the table
        addDictToTab({"Chromosome1": record.CHROM if "chr" not in record.CHROM else record.CHROM[3:],
                     "Position1": record.POS, "ID": record.ID, "Reference": record.REF, "Alternate": record.ALT[0],
                     "Quality": record.QUAL if record.QUAL is not None else ".", 
                     "Filter": "PASS" if not "LowQual" in record.FILTER else "LowQual",
                     "Sketchy?": sketchy})
        
        # Processes the format information in the .vcf into a dictionary,
        # for each sample
        fields = record.FORMAT.split(":")
        
        # If two .bam files are given, assume they're a normal and a tumor.
        # Otherwise just label them as "samples"
        if(len(record.samples) == 2):
                nKeys = ["NORMAL" + "_" + i for i in fields]
                nVals = [record.samples[0][field] for field in fields]
                newValues = dict(zip(nKeys, nVals))
                addDictToTab(newValues)

                tKeys = ["TUMOR" + "_" + i for i in fields]
                tVals = [record.samples[1][field] for field in fields]
                newValues = dict(zip(tKeys, tVals))
                addDictToTab(newValues)
        else:
            for sample in range(0, len(record.samples)):
                sKeys = ["SAMPLE" + str(sample + 1) + "_" + i for i in fields]
                sVals = [record.samples[sample][field] for field in fields]
                newValues = dict(zip(sKeys, sVals))
            
                # Add the new sample information to the .tab file
                addDictToTab(newValues)     
        
        # Add in the name of the genes at the breakpoints using the
        # geneName module, as well as the genes in between
        bps = {"Breakpoint1" : geneName(row[tabFile[0].index("Chromosome1")], 
                                        int(row[tabFile[0].index("Position1")]),
                                        ucsc_info_path)}
        # If its in the file, find the second gene breakpoint too
        if(row[tabFile[0].index("Chromosome2")] != NULL_VAL):
               bps["Breakpoint2"] = geneName(row[tabFile[0].index("Chromosome2")],
                                         int(row[tabFile[0].index("Position2")]),
                                         ucsc_info_path)
               # If its on the same gene, find the genes between too
               if row[tabFile[0].index("Chromosome1")] == row[tabFile[0].index("Chromosome2")]:
                       bps["Between"] = geneBetween(row[tabFile[0].index("Chromosome2")],
                                                    int(row[tabFile[0].index("Position1")]),
                                                    int(row[tabFile[0].index("Position2")]),
                                                    MAX_LIST,
                                                    ucsc_info_path)
                       bps["Exons"] = exonBetween(row[tabFile[0].index("Chromosome2")],
                                            int(row[tabFile[0].index("Position1")]),
                                            int(row[tabFile[0].index("Position2")]),
                                            MAX_LIST,
                                            ucsc_info_path)
                       # bps["Common"] = polymorphisms(row[tabFile[0].index("Chromosome2")],
                       #                      int(row[tabFile[0].index("Position1")]),
                       #                      int(row[tabFile[0].index("Position2")]))
                       bps["Common"] = np.nan
               else:
                       bps["Exons"] = (exonName(row[tabFile[0].index("Chromosome1")],
                                                int(row[tabFile[0].index("Position1")]),
                                                ucsc_info_path) +
                                       exonName(row[tabFile[0].index("Chromosome2")],
                                                int(row[tabFile[0].index("Position2")]),
                                                ucsc_info_path))
        
        # Calculate svlen and add it to the .tab file
        svlen = {"SVLEN" : (record.INFO["END"] - record.POS)}
        addDictToTab(svlen)

        # Add the breakpoint names to the .tab file
        addDictToTab(bps)

        # Add the newly aligned row into the final array
        tabFile.append(row)

    # Fill in the empty spots with NULL_VAL, so that the columns align properly
    table_length = len(tabFile[0])
    for row in tabFile:
        while len(row) < table_length:
            row.append(NULL_VAL)

    # Write the data to a new file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(tabFile)

    if (filter_info_file is not None) and (sv_type is not None):
        with open(filter_info_file, 'w') as filter_info:
            filter_info.write(sv_type + "\tunfiltered\t%d\n" % unfiltered_count)
            filter_info.write(sv_type + "\tvcf2tab_qualilty_filters\t%d\n" % qual_filtered_count)

    print("Finished converting .vcf file to .tab file.")

if __name__ == '__main__':

    # Process the command line input using optparse
    parser = OptionParser()
    parser.add_option("--vcf", help="the location of the .vcf file to process", dest="vcf")    
    parser.add_option("--bam", help="the path to the bam file, used for marking sketchy regions", dest="bam", default=".")
    parser.add_option("--lowqual", help="the quality of filtering to apply to the vcf file", dest="lowqual", default=False)
    parser.add_option("--output_file", help="The file to store output in", dest="output_file")
    parser.add_option("--ucsc_info_path", help="The directory containing the UCSC_info directory", dest="ucsc_info_path")
    parser.add_option("--filter_info_file", help = "The file to save filtering info to", dest="filter_info_file")
    parser.add_option("--sv_type", help = "The SV type", dest = "sv_type")
    (options, args) = parser.parse_args()

    options.lowqual = bool(options.lowqual)

    # Raise optparse errors as necessary
    if options.vcf is None or options.vcf[-4:] != ".vcf":
        parser.error("An invalid .vcf file was given")
    if not(os.path.isfile(options.vcf)):
        parser.error("The .vcf file given does not exist") 
    elif not options.bam == "." and options.bam[-4:] != ".bam":
        parser.error("An invalid .bam file was given")
    if options.lowqual != True and options.lowqual != False:
        parser.error("An invalid lowqual value was given")
    if options.ucsc_info_path is None:
        parser.error("Please provide ucsc_info_path")
    if options.output_file is None:
        parser.error("Please provide output_file")

    # If no errors were raised, let the processing begin
    vcf2Tab(options.vcf, options.lowqual, options.bam, options.output_file, options.ucsc_info_path, options.filter_info_file, options.sv_type)
