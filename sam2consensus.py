#!/usr/bin/env python2
# -*- coding: utf-8 -*-



'''
The program takes as input a SAM file resulting from mapping short reads to a collection of
gene sequences as reference, then it calculates the consensus sequence per gene without
considering the reference. It adds insertions and can take a custom consensus threshold,
the consensus method is the same as the one described for Geneious
(http://assets.geneious.com/manual/8.1/GeneiousManualse41.html). Regions with no coverage
are filled with Ns.

Input SAM files have to be sorted, and contain only mapped reads (preferably, to save
computation times). Original reference FASTAs are not necessary.

It will produce a FASTA sequence per gene, and in case the gene has insertions it will also create
a separate SAM file just for the particular gene for verification purposes.
'''


__author__      = "Edgardo M. Ortiz"
__credits__     = "Deise J.P. Gonçalves"
__version__     = "1.5"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2017-10-07"


import sys
import re
import operator
import argparse
import gzip
import io

def parsecigar(cigarstring, seq, pos_ref):
    '''
    Modifies a sequence according to its CIGAR string

    :param cigarstring: cigar string, check SAM format specification
    :param seq: raw sequence
    :param pos_ref: position in the reference of the leftmost aligned nucleotide
    :return: edited sequence according to cigar string, list of tuples for
             insertions indicating coordinate in the reference and the
             sequence inserted
    '''

    matches = re.findall(r"(\d+)([MIDNSHPX=]{1})", cigarstring)
    cigar = [{"type": m[1], "length": int(m[0])} for m in matches]
    start = 0
    start_ref = pos_ref
    seqout = ""
    insert = []
    for c in range(0, len(cigar)):
        l = cigar[c]["length"]
        if cigar[c]["type"] == "S":
            start += l
        elif cigar[c]["type"] == "H":
            continue
        elif cigar[c]["type"] == "M":
            seqout += seq[start:start + l]
            start += l
            start_ref += l
        elif cigar[c]["type"] == "=":
            seqout += seq[start:start + l]
            start += l
            start_ref += l
        elif cigar[c]["type"] == "X":
            seqout += seq[start:start + l]
            start += l
            start_ref += l
        elif cigar[c]["type"] == "I":
            insert.append((start_ref, seq[start:start+l]))
            start += l
        elif cigar[c]["type"] == "D":
            seqout += "-" * l
            start_ref += l
        elif cigar[c]["type"] == "N":
            seqout += "-" * l
            start_ref += l
        elif cigar[c]["type"] == "P": # check functionality of this case
            seqout += "-" * l
            start_ref += l
        else:
            print "SAM file probably contains unmapped reads"
    return seqout, insert


def main():
    parser = argparse.ArgumentParser(description="Calculates the consensus sequence from reads aligned to a multi-gene fasta reference")
    parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
        help="Name of the SAM file, SAM must be sorted and can be compressed with gzip")
    parser.add_argument("-c", "--consensus-threshold", action="store", dest="cons_threshold", type=float, default=0.25,
        help="Consensus threshold sensu Geneious, default=0.25")
    parser.add_argument("-m", "--min-depth", action="store", dest="min_depth", type=int, default=5,
    	help="Minimum read depth at each site to report the nucleotide in the consensus, default=5")
    parser.add_argument("-o", "--outfolder", action="store", dest="outfolder", default="./",
        help="Name of output folder, default=same folder as input")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="",
    	help="Prefix for output file name, default=input filename without .sam extension")
    parser.add_argument("-f", "--fill", action="store", dest="fill", default="-",
        help="Character for padding regions not covered in the reference, default= - (gap)")
    parser.add_argument("-s", "--sam-verify", action="store_true", dest="sam_verify", default=False,
        help="Enable creation of individual SAM files per gene for verification of insertion, disabled by default")
    args = parser.parse_args()

    filename = args.filename
    cons_threshold = args.cons_threshold
    if args.prefix == "":
    	prefix = ".".join(args.filename.split(".")[:-1]).split("/")[-1]   # Input filename without extension
    else:
    	prefix = args.prefix
    outfolder = args.outfolder

    min_depth = args.min_depth

    fill = args.fill

    sam_verify = args.sam_verify

    print "Processing sample "+prefix+" ...\n"
    
    if outfolder[-1] != "/":
        outfolder += "/"
    
    if filename.endswith(".gz"):
        opener = gzip.open(filename, "rb")
    else:
        opener = open(filename, "r")

    # Process the SAM file in a single pass
    with opener as mapfile:
        genes = {}                                              # Container of sequences per gene
        insertions = {}                                         # Container for insertions with coordinates per gene
        gene_previous = ""                                      # Stores name of previous gene processed
        sam_file = []                                           # Container for SAM file whe insertions are present
        sam_header = ""                                         # Store the SAM header
        sam_reads = []                                          # Container for reads to be written to the SAM

        for line in mapfile:
        	# Store SAM header
            if line.startswith("@HD"):
                sam_header = line.strip("\n")
                sam_file.append(line.strip("\n"))

            # Extract gene names
            elif line.startswith("@SQ"):

                # Obtain the name of the first gene in the file
                gene_name = line.split("\t")[1].replace("SN:","").split()[0] # Get gene name, skip description if present

                if gene_previous == "":
                    gene_previous = gene_name

                # Populate empty dictionary, values to be stored in a list per gene
                genes[gene_name] = []
                insertions[gene_name] = []

                # Populate each gene with as many empty nucleotides as the reference
                for nuc in range(0, int(line.split("\t")[2].replace("LN:",""))):
                    genes[gene_name].append({"A":0,"C":0,"T":0,"G":0,"-":0,"N":0})

                # Also add the length of each gene after the list of nucleotides
                genes[gene_name].append(int(line.split("\t")[2].replace("LN:","")))

            # Start processing the aligned reads, skip unaligned [*]
            elif line[0] != "@" and line.split("\t")[5] != "*":
                gene_current = line.split("\t")[2].split()[0] # Get gene name, skip description if present

                # If we haven't started processing the next gene...
                if gene_current == gene_previous:
                    pos_ref = int(line.split("\t")[3]) - 1              # Starting position in the reference
                    cigar = line.split("\t")[5]                         # CIGAR string fo the aligned read
                    seqraw = line.split("\t")[9]                        # Unaltered sequence of the aligned read
                    seqout, insert = parsecigar(cigar, seqraw, pos_ref) # Parse the CIGAR and obtain edited sequence
                                                                        # and list of insertions
                    
                    # Fill the nucleotides in the respective gene according to the
                    # sequence processed according to its CIGAR string
                    for pos in range(0,len(seqout)):
                        genes[gene_current][pos+pos_ref][seqout[pos]] += 1

                    # Add insertions with coordinates to the dictionary of insertions per gene
                    for ins in insert:
                        insertions[gene_current].append(ins)

                    # Update name of previous gene with current
                    gene_previous = gene_current

                    # Add the read in case a SAM is produced for this gene
                    sam_reads.append(line.strip("\n"))

                # If we started processing the next gene, stop and summarize the previous gene,
                # then continue as normal for next gene...
                else:
                    # Calculate average coverage per base and add to the gene info
                    nuc_covs = 0
                    for pos in range(0, (len(genes[gene_previous])-1)):
                        nuc_covs += sum(genes[gene_previous][pos].values())
                    cov_average = float(nuc_covs)/float(genes[gene_previous][-1])
                    genes[gene_previous].append(cov_average)

                    # Find real insertions based on coverage of adjacent nucleotides
                    real_insertions_coordinates = []
                    real_insertions_motifs = []
                    for ins in sorted(set(insertions[gene_previous])):

                        # Get the average coverage of the nucleotide before and after the insertion
                        if genes[gene_previous][-2]-1 > ins[0]:
                            cov_at_edges = float(sum(genes[gene_previous][ins[0]].values())+sum(genes[gene_previous][ins[0]+1].values()))/2
                        else:
                            cov_at_edges = float(sum(genes[gene_previous][ins[0]].values()))

                        # If the insertion has acceptable coverage accept it as real
                        if cov_at_edges >= min_depth:
	                        if (insertions[gene_previous].count(ins) >= (cov_at_edges-insertions[gene_previous].count(ins))) and (insertions[gene_previous].count(ins) >= cov_at_edges*cons_threshold):
	                            real_insertions_coordinates.append(ins[0])
	                            real_insertions_motifs.append(ins[1])
	                            print "Insertion detected in gene "+gene_previous+", coverage at sides of insertion: "+str(cov_at_edges)+", coverage of the insertion: "+str(insertions[gene_previous].count(ins))+", position/motif: "+str(ins)

                    # If the gene has real insertions produce a SAM for verification and
                    # eliminate insertions with low coverage (errors)
                    if real_insertions_coordinates != []:
                        if sam_verify == True:
                            print gene_previous+" contains insertion(s), a separate SAM file will be additionally created for this gene."
                            sam_file.append("@SQ"+"\t"+"SN:"+gene_previous+"\t"+"LN:"+str(genes[gene_previous][-2]))
                            for read in sam_reads:
                                sam_file.append(read)
                            outfile = open(outfolder+gene_previous+"_to_"+prefix+".sam", "w")
                            outfile.write("\n".join(sam_file)+"\n")
                        insertions[gene_previous] = [real_insertions_coordinates,real_insertions_motifs]
                    else:
                        del insertions[gene_previous]
                    
                    # Reset SAM header, empty list of SAM reads
                    sam_file = [sam_header]
                    sam_reads = []
                    print "Gene "+gene_previous+" processed\n"

                    # Process current read as in line 92
                    pos_ref = int(line.split("\t")[3]) - 1
                    cigar = line.split("\t")[5]
                    seqraw = line.split("\t")[9]
                    seqout, insert = parsecigar(cigar, seqraw, pos_ref)
                    for site in range(0,len(seqout)):
                        genes[gene_current][site+pos_ref][seqout[site]] += 1
                    for ins in insert:
                        insertions[gene_current].append(ins)
                    gene_previous = gene_current
                    sam_reads.append(line.strip("\n"))

        # For the last read of the last gene only, 
        # same process as line 117
        nuc_covs = 0
        for pos in range(0, (len(genes[gene_current])-1)):
            nuc_covs += sum(genes[gene_current][pos].values())
        cov_average = float(nuc_covs)/float(genes[gene_current][-1])
        genes[gene_current].append(cov_average)

        real_insertions_coordinates = []
        real_insertions_motifs = []
        for ins in sorted(set(insertions[gene_current])):
            if genes[gene_current][-2]-1 > ins[0]:
                cov_at_edges = float(sum(genes[gene_current][ins[0]].values())+sum(genes[gene_current][ins[0]+1].values()))/2
            else:
                cov_at_edges = float(sum(genes[gene_current][ins[0]].values()))

            if cov_at_edges >= min_depth:
	            if (insertions[gene_current].count(ins) >= (cov_at_edges-insertions[gene_current].count(ins))) and (insertions[gene_current].count(ins) >= cov_at_edges*cons_threshold):
	                real_insertions_coordinates.append(ins[0])
	                real_insertions_motifs.append(ins[1])
	                print "Insertion detected in gene "+gene_current+", coverage at sides of insertion: "+str(cov_at_edges)+", coverage of the insertion: "+str(insertions[gene_current].count(ins))+", position/motif: "+str(ins)
        
        if real_insertions_coordinates != []:
            if sam_verify == True:
                print gene_current+" contains insertion(s), a separate SAM file will be additionally created for this gene."
                sam_file.append("@SQ"+"\t"+"SN:"+gene_current+"\t"+"LN:"+str(genes[gene_current][-2]))
                for read in sam_reads:
                    sam_file.append(read)
                outfile = open(outfolder+gene_current+"_to_"+prefix+".sam", "w")
                outfile.write("\n".join(sam_file)+"\n")
            insertions[gene_current] = [real_insertions_coordinates,real_insertions_motifs]
        else:
            del insertions[gene_current]

        print "Gene "+gene_current+" processed\n"



    # Dictionary to translate IUPAC ambiguities
    amb = {("-","A"):"A",
           ("-","C"):"C",
           ("-","G"):"G",
           ("-","N"):"-",
           ("-","T"):"T",
           ("A","C"):"M",
           ("A","G"):"R",
           ("A","N"):"A",
           ("A","T"):"W",
           ("C","G"):"S",
           ("C","N"):"C",
           ("C","T"):"Y",
           ("G","N"):"G",
           ("G","T"):"K",
           ("N","T"):"T"}



    # If a gene had zero reads mapped then set its coverage to 0.0 
    for gene in genes:
        if type(genes[gene][-2]) != int:
            genes[gene].append(float(0))



    # Obtain sequence from the "genes" dictionary
    fastas = {}
    for gene in genes:
        for pos in range(0, genes[gene][-2]):
            count_nucs = list(sorted(genes[gene][pos].iteritems(), key=operator.itemgetter(1), reverse=True)[:2])
            cov_site = float(sum(genes[gene][pos].values()))
            if cov_site < min_depth:
                if gene not in fastas:
                    fastas[gene] = fill
                else:
                    fastas[gene] += fill
            elif count_nucs[0][1] >= cons_threshold*cov_site:
                if gene not in fastas:
                    fastas[gene] = count_nucs[0][0]
                else:
                    fastas[gene] += count_nucs[0][0]
            else:
                if gene not in fastas:
                    fastas[gene] = amb[tuple(sorted((count_nucs[0][0],count_nucs[1][0])))]
                else:
                    fastas[gene] += amb[tuple(sorted((count_nucs[0][0],count_nucs[1][0])))]



    # Add insertions for genes with real insertions
    if insertions:
        for gene in insertions:
            if insertions[gene] != []:
                seq = fastas[gene]
                start = 0
                seqout = ""
                for i in range(0, len(insertions[gene][0])):
                    seqout += seq[start:insertions[gene][0][i]]
                    seqout += insertions[gene][1][i]
                    start = insertions[gene][0][i]
                seqout += seq[insertions[gene][0][-1]:]
                fastas[gene] = seqout



    # Write fasta output files
    for gene in fastas:
        outfile = open(outfolder+prefix+"_to_"+gene+"_cons"+str(cons_threshold)+".fasta", "w")
        outfile.write(">"+prefix+" Mapped to: "+gene+", consensus threshold: "+str(cons_threshold)+", coverage: "+str(genes[gene][-1])+"\n"+fastas[gene]+"\n")



if __name__ == "__main__":
    main()
