#! python3

# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.0                 ?                    ?        initial development
#     0.1        2018-12-17        Drew Thompson            updated for WDL
#     0.2        2019-07-31        Drew Thompson      remove harcoded paths

# UCSC_naming is a module written to find the name of a gene in a given
# interval. geneName takes a chromosome number, a position, and an
# uncertainty in the position, and returns any gene that happens to be
# in that interval.

import os
import sys
import csv
 
GENE_NAMES_BED = "GENE_NAMES"     # Path to gene .bed file
EXON_NAMES_BED = "EXON_NAMES" # Path to exon .bed file
COMMON_POLYMORPHISMS_BED = "COMMON_POLYMORPHISMS"
NULL_VAL = "NONE"                 # Value to return if no gene found
EXON_GIVE = 50
POLY_UNCERTAINTY = 100

def polymorphisms(chromosome, startI, endI, ucsc_info_path):
    """
    Returns a list of polymorphisms that appear in the specified 
    interval for healthy people

    Args:
        chromosome (string): The chromosome to search on
        start (int): The start point of the interval
        end (int): The end point of the interval

    Returns:
        A list of strings, each of which corresponds to a common polymorphism
    """

    # Load the mutations and prepare for binary searching
    PACKAGE_DIR = "%s/UCSC_info/%s.bed" % (ucsc_info_path,
                                    COMMON_POLYMORPHISMS_BED)
    polymorphs = list(csv.reader(open(PACKAGE_DIR, 'rt'), delimiter='\t'))
    start, end, mid, point = 0, len(polymorphs), 0, startI
    prev = -1

    # Start the binary search
    while True:
        if abs(start - end) <= 2:
            break
        mid = int((start + end) / 2)
        if mid == prev:
            break
        if point > int(polymorphs[mid][2]):
            start = mid
        elif point < int(polymorphs[mid][1]):
            end = mid
        prev = mid

    intervals, up, down = [], mid, mid + 1

    # Process all upstream 
    while -POLY_UNCERTAINTY <= int(polymorphs[up][1]) - startI <= POLY_UNCERTAINTY:
         if (polymorphs[up][0] == str(chromosome) and
             -POLY_UNCERTAINTY <= int(polymorphs[up][2]) - endI <= POLY_UNCERTAINTY):
             intervals.append(up)
         up -= 1

    # Process all downstream
    while -POLY_UNCERTAINTY <= int(polymorphs[down][1]) - startI <= POLY_UNCERTAINTY:
         if (polymorphs[down][0] == str(chromosome) and
             -POLY_UNCERTAINTY <= int(polymorphs[down][2]) - endI <= POLY_UNCERTAINTY):
             intervals.append(down)
         down += 1

    # Process the findings
    intervals = [polymorphs[i][3] for i in intervals]
    return list(set(intervals))

def exonBetween(chromosome, point1, point2, max, ucsc_info_path):
    """
    exonBetween(chromosome, point) takes in a chromsome number and
    start and end, and returns the names of the genes at that location.

    Args:
	chromosome (string): The chromosome number being searched for
	point1 (int): The approximate start location of the exon
    point2 (int): The approximate end location of the exon
    max (int): The maximum amount of genes to return before shortening

    Returns:
    	string: The name of the exons at that location. 
    """ 

    # Load the genes and prepare for binary searching
    if "chr" in chromosome:
        chromosome = chromosome[3:]
    PACKAGE_DIR = "%s/UCSC_info/%s_%s.bed" % (ucsc_info_path,
                                   EXON_NAMES_BED, chromosome)
    exons = list(csv.reader(open(PACKAGE_DIR, 'rt'), delimiter='\t'))
    start, end, mid, prev = 0, len(exons), 0, -1
    point = (point1 + point2) / 2

    # Checks if the intervals intersect
    def captured(exon):
        if ((point1 <= int(exons[exon][2]) <= point2) or
            (point1 <= int(exons[exon][1]) <= point2) or
            (int(exons[exon][1]) < point1 and int(exons[exon][2]) > point2) or
            (int(exons[exon][1]) > point1 and int(exons[exon][2]) < point2)):
            return True
        else:
            return False

    # Start the binary search
    while True:
        mid = int((start + end) / 2)
        if mid == prev:
            break
        elif point > int(exons[mid][2]):
            start = mid
        elif point < int(exons[mid][1]):
            end = mid
        prev = mid
     
    # Find all genes in this area
    names = ([mid] if captured(mid) else [])
    downstream, upstream = mid - 1, mid + 1
    while downstream > -1 and captured(downstream): 
        names.append(downstream)
        downstream -= 1
    names = names[::-1]
    while upstream < len(exons) and captured(upstream):
         names.append(upstream)
         upstream += 1
    names = [exons[i] for i in names]

    # Report back the gene found
    if len(names) > max:
        return ["%s-%s" % (names[0][3], names[0][4]), "...", "%s%s" % (names[-1][3], names[-1][4])]
    else:
        return list(set(["%s-%s" % (exon[3], exon[4]) for exon in names]))

def geneBetween(chromosome, point1, point2, max, ucsc_info_path):
    """
    geneBetween(chromosome, point) takes in a chromsome number and
    start and end, and returns the names of the genes at that location.

    Args:
	chromosome (string): The chromosome number being searched for
	point1 (int): The approximate start location of the gene
    point2 (int): The approximate end location of the gene
    max (int): The maximum amount of genes to return before shortening

    Returns:
    	string: The name of the gene at that location. 
    """ 

    # Load the genes and prepare for binary searching
    if "chr" in chromosome:
        chromosome = chromosome[3:]
    PACKAGE_DIR = "%s/UCSC_info/%s_%s.bed" % (ucsc_info_path,
                                   GENE_NAMES_BED, chromosome)
    genes = list(csv.reader(open(PACKAGE_DIR, 'rt'), delimiter='\t'))
    start, end, mid, prev = 0, len(genes), 0, -1
    point = (point1 + point2) / 2

    # Checks if the intervals intersect
    def captured(gene):
        if ((point1 <= int(genes[gene][2]) <= point2) or
            (point1 <= int(genes[gene][1]) <= point2) or
            (int(genes[gene][1]) < point1 and int(genes[gene][2]) > point2) or
            (int(genes[gene][1]) > point1 and int(genes[gene][2]) < point2)):
            return True
        else:
            return False

    # Start the binary search
    while True:
        mid = int((start + end) / 2)
        if mid == prev:
            break
        elif point > int(genes[mid][2]):
            start = mid
        elif point < int(genes[mid][1]):
            end = mid
        prev = mid
     
    # Find all genes in this area
    names = ([mid] if captured(mid) else [])
    downstream, upstream = mid - 1, mid + 1
    while downstream > -1 and captured(downstream): 
        names.append(downstream)
        downstream -= 1
    names = names[::-1]
    while upstream < len(genes) and captured(upstream):
         names.append(upstream)
         upstream += 1
    names = [genes[i] for i in names]

    # Report back the gene found
    if len(names) > max:
        return ["%s%s" % (names[0][3], names[0][4]), "...", "%s%s" % (names[-1][3], names[-1][4])]
    else:
        return ["%s%s" % (gene[3], gene[4]) for gene in names]

def exonName(chromosome, point, ucsc_info_path):
    """
    exonName(chromosome, point) takes in a chromsome number and
    breakpoint location, and returns the name of the exon at that location.


    chromosome (string): The chromosome number being searched for
    point (int): The approximate location of the exon

    Returns:
        string: The name of the gene at that location. 
    """

    # Load the genes and prepare for binary searching
    if "chr" in chromosome:
        chromosome = chromosome[3:]
    PACKAGE_DIR = "%s/UCSC_info/%s_%s.bed" % (ucsc_info_path,
                                   EXON_NAMES_BED, chromosome)
    exons = list(csv.reader(open(PACKAGE_DIR, 'rt'), delimiter='\t'))
    start, end, mid, prev = 0, len(exons), 0, -1

    # Start the binary search
    while True:
        mid = int((start + end) / 2)
        if mid == prev:
            break
        if point > int(exons[mid][2]):
            start = mid
        elif point < int(exons[mid][1]):
            end = mid
        prev = mid

    # Find all genes in this area
    names = ([mid] if (int(exons[mid][1]) <= point <= int(exons[mid][2])) else [])
    downstream, upstream = mid - 1, mid + 1
    while downstream < len(exons) and (int(exons[downstream][1]) <= point <= int(exons[downstream][2])):
        names.append(downstream)
        downstream -= 1
    names = names[::-1]
    while upstream < len(exons) and (int(exons[upstream][1]) <= point <= int(exons[upstream][2])):
         names.append(upstream)
         upstream += 1
    names = [exons[i] for i in names]
    if names == []:
        if point <= int(exons[mid][1]):
            upstream = mid
        else:
            downstream = mid
    if upstream > len(exons) - 1:
        upstream = len(exons)
        exons.append(["NULL", point, point, "RIGHT_END", 0, "+"])
    if downstream < 0:
        downstream = 0
        exons = [["NULL", point, point, "LEFT_END", 0, "+"]] + exons

    # Report back the gene found
    if len(names) > 0:
        return list(set(["%s-%s" % (exon[3], exon[4]) for exon in names]))

    # If none were found, get the nearest upstream and downstream genes (look
    # until the breakpoints don't match anymore, and then stop, starting from
    # the location previously stopped at).
    else:
        names = [downstream, upstream]
        return (["(%s-%s,%s)" % (exons[exon][3], exons[exon][4], min(abs(int(exons[exon][1]) - point),
                                                                     abs(int(exons[exon][2]) - point)))
                for exon in names])


def geneName(chromosome, point, ucsc_info_path):
    """
    geneName(chromosome, point) takes in a chromsome number and
    breakpoint location, and returns the name of the gene at that location.

    Args:
	chromosome (string): The chromosome number being searched for
	point (int): The approximate start location of the gene

    Returns:
    	string: The name of the gene at that location. 
    """ 

    # Load the genes and prepare for binary searching
    if "chr" in chromosome:
        chromosome = chromosome[3:]
    PACKAGE_DIR = "%s/UCSC_info/%s_%s.bed" % (ucsc_info_path,
                                   GENE_NAMES_BED, chromosome)
    genes = list(csv.reader(open(PACKAGE_DIR, 'rt'), delimiter='\t'))
    start, end, mid, prev = 0, len(genes), 0, -1

    # Start the binary search
    while True:
        mid = int((start + end) / 2)
        if mid == prev:
            break
        if point > int(genes[mid][2]):
            start = mid
        elif point < int(genes[mid][1]):
            end = mid
        prev = mid

    # Find all genes in this area
    names = ([mid] if (int(genes[mid][1]) <= point <= int(genes[mid][2])) else [])
    downstream, upstream = mid - 1, mid + 1
    while downstream < len(genes) and (int(genes[downstream][1]) <= point <= int(genes[downstream][2])):
        names.append(downstream)
        downstream -= 1
    names = names[::-1]
    while upstream < len(genes) and (int(genes[upstream][1]) <= point <= int(genes[upstream][2])):
         names.append(upstream)
         upstream += 1
    names = [genes[i] for i in names]
    if names == []:
        if point <= int(genes[mid][1]):
            upstream = mid
        else:
            downstream = mid 
    if upstream > len(genes) - 1:
        genes = genes + [["NULL", point, point, "RIGHT_END", 0, "+"]]
        upstream = len(genes) - 1
    if downstream < 0:
        genes = [["NULL", point, point, "LEFT_END", 0, "+"]] + genes
        downstream = 0

    # Report back the gene found
    if len(names) > 0:
        return list(set(["%s%s" % (gene[3], gene[4]) for gene in names]))

    # If none were found, get the nearest upstream and downstream genes (look
    # until the breakpoints don't match anymore, and then stop, starting from
    # the location previously stopped at).
    else:
        names = [downstream, upstream]
        return (["(%s%s,%s)" % (genes[gene][3], genes[gene][4], min(abs(int(genes[gene][1]) - point), 
                                                                    abs(int(genes[gene][2]) - point))) 
                for gene in names])
