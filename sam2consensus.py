#!/usr/bin/env python2
# -*- coding: utf-8 -*-



'''
The program takes as input a SAM file resulting from mapping short reads to a collection of contigs as reference (the
contigs can correspond to separate genes for example), then it calculates the consensus sequence per contig without
considering the reference. It adds insertions and can take a custom consensus threshold, the consensus method is the
same as the one described for Geneious (http://assets.geneiious.com/manual/8.1/GeneiousManualse41.html).
Regions with no coverage are filled with Ns.

Input SAM files have to be sorted. Original reference FASTAs are not necessary since the consensus is reference-free.

It will produce a FASTA sequence per contig, and in case the contig has insertions it can also create
a separate SAM file just for the particular contig for verification purposes.
'''



__author__      = "Edgardo M. Ortiz"
__credits__     = "Deise J.P. Gonçalves"
__version__     = "1.6"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2018-08-09"



import sys
import re
import operator
import argparse
import gzip



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
	parser = argparse.ArgumentParser(description="Calculates the consensus sequence from reads aligned to a multi-contig fasta reference")
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
	parser.add_argument("-n", action="store", dest="n", default=80,
		help="Split FASTA output sequences every n nucleotides, default=80")
	parser.add_argument("-s", "--sam-verify", action="store_true", dest="sam_verify", default=False,
		help="Enable creation of individual SAM files per contig for verification of insertion, disabled by default")
	args = parser.parse_args()

	filename = args.filename
	# Prepare the opener if the SAM is compressed
	if filename.endswith(".gz"):
		opener = gzip.open
	else:
		opener = open

	cons_threshold = args.cons_threshold
	
	if args.prefix == "":
		prefix = "".join(args.filename.split("/")[-1]).split(".")[0]   # Prefix will be input filename without extension
	else:
		prefix = args.prefix
	
	outfolder = args.outfolder
	if outfolder[-1] != "/":
		outfolder += "/"

	min_depth = args.min_depth

	nchar = args.n

	fill = args.fill

	sam_verify = args.sam_verify



	# Start processing SAM file
	print "Processing sample "+prefix+"...\n"



	# First process header of file to extract info about contigs
	with opener(filename) as mapfile:
		contigs = {}                                            # Container of sequences per contig
		insertions = {}                                         # Container for insertions with coordinates per contig
		contig_previous = ""                                    # Stores name of previous contig processed
		sam_file = []                                           # Container for SAM file when insertions are present
		sam_header = ""                                         # Store the SAM header
		sam_reads = []                                          # Container for reads to be written to the SAM

		for line in mapfile:
			# Store SAM header
			if line.startswith("@HD"):
				sam_header = line.strip("\n")
				sam_file.append(line.strip("\n"))

			# Extract contig names
			elif line.startswith("@SQ"):

				# Obtain the name of the first contig in the file
				contig_name = line.split("\t")[1].replace("SN:","").split()[0] # Get contig name, skip description if present

				if contig_previous == "":
					contig_previous = contig_name

				# Populate empty dictionary, values to be stored in a list per contig
				contigs[contig_name] = []
				insertions[contig_name] = []

				# Populate each contig with as many empty nucleotides as the reference
				for nuc in range(0, int(line.split("\t")[2].replace("LN:",""))):
					contigs[contig_name].append({"A":0,"C":0,"T":0,"G":0,"-":0,"N":0})

				# Also add the length of each contig after the list of nucleotides
				contigs[contig_name].append(int(line.split("\t")[2].replace("LN:","")))

			elif not line.startswith("@"):
				break

	mapfile.close()



	# Now process the reads that were aligned to the reference contigs
	with opener(filename) as mapfile:
		# Set  read counter to show progress
		read_count = 0

		while True:
			# Load big file by chunks;
			mapfile_chunk = mapfile.readlines(100000)
			if not mapfile_chunk:
				break

			for line in mapfile_chunk:
				if line[0] != "@" and line.split("\t")[5] != "*":
					read_count += 1
					line_split = line.split("\t")

					if read_count % 1000000 == 0:
						print str(read_count)+" reads processed"

					contig_current = line_split[2].split()[0] # Get contig name, skip description if present

					# If we haven't started processing the next contig...
					if contig_current == contig_previous:
						
						pos_ref = int(line_split[3]) - 1              # Starting position in the reference
						# cigar  = line_split[5]                      # CIGAR string fo the aligned read
						# seqraw = line_split[9]                      # Unaltered sequence of the aligned read
						# Parse the CIGAR and obtain edited sequence and list of insertions
						seqout, insert = parsecigar(line_split[5], line_split[9], pos_ref)
						
						# Fill the nucleotides in the respective contig according to the
						# sequence processed according to its CIGAR string
						for pos in range(0,len(seqout)):
							contigs[contig_current][pos+pos_ref][seqout[pos]] += 1

						# Add insertions with coordinates to the dictionary of insertions per contig
						for ins in insert:
							insertions[contig_current].append(ins)

						# Update name of previous contig with current
						contig_previous = contig_current

						# Add the read in case a SAM is produced for this contig
						if sam_verify == True:
							sam_reads.append(line.strip("\n"))

					# If we started processing the next contig, stop and summarize the previous contig,
					# then continue as normal for next contig...
					else:
						# Calculate average coverage per base and add to the contig info
						nuc_covs = 0
						for pos in range(0, (len(contigs[contig_previous])-1)):
							nuc_covs += sum(contigs[contig_previous][pos].values())
						cov_average = float(nuc_covs)/float(contigs[contig_previous][-1])
						contigs[contig_previous].append(cov_average)

						# Find real insertions based on coverage of adjacent nucleotides
						real_insertions_coordinates = []
						real_insertions_motifs = []
						for ins in sorted(set(insertions[contig_previous])):

							# Get the average coverage of the nucleotide before and after the insertion
							if contigs[contig_previous][-2]-1 > ins[0]:
								cov_at_edges = float(sum(contigs[contig_previous][ins[0]].values())+sum(contigs[contig_previous][ins[0]+1].values()))/2
							else:
								cov_at_edges = float(sum(contigs[contig_previous][ins[0]].values()))

							# If the insertion has acceptable coverage accept it as real
							if cov_at_edges >= min_depth:
								if (insertions[contig_previous].count(ins) >= (cov_at_edges-insertions[contig_previous].count(ins))) and (insertions[contig_previous].count(ins) >= cov_at_edges*cons_threshold):
									real_insertions_coordinates.append(ins[0])
									real_insertions_motifs.append(ins[1])
									print "Insertion detected in contig "+contig_previous+", coverage at sides of insertion: "+str(cov_at_edges)+", coverage of the insertion: "+str(insertions[contig_previous].count(ins))+", position/motif: "+str(ins)

						# If the contig has real insertions produce a SAM for verification and
						# eliminate insertions with low coverage (errors)
						if real_insertions_coordinates != []:
							if sam_verify == True:
								print contig_previous+" contains insertion(s), a separate SAM file will be additionally created for this contig."
								sam_file.append("@SQ"+"\t"+"SN:"+contig_previous+"\t"+"LN:"+str(contigs[contig_previous][-2]))
								for read in sam_reads:
									sam_file.append(read)
								outfile = open(outfolder+contig_previous+"_to_"+prefix+".sam", "w")
								outfile.write("\n".join(sam_file)+"\n")
							insertions[contig_previous] = [real_insertions_coordinates,real_insertions_motifs]
						else:
							del insertions[contig_previous]
						
						# Reset SAM header, empty list of SAM reads
						sam_file = [sam_header]
						sam_reads = []
						print "Contig "+contig_previous+" processed\n"

						# Continue normal processing of reads
						pos_ref = int(line.split("\t")[3]) - 1
						cigar = line.split("\t")[5]
						seqraw = line.split("\t")[9]
						seqout, insert = parsecigar(cigar, seqraw, pos_ref)
						for site in range(0,len(seqout)):
							contigs[contig_current][site+pos_ref][seqout[site]] += 1
						for ins in insert:
							insertions[contig_current].append(ins)
						contig_previous = contig_current
						if sam_verify == True:
							sam_reads.append(line.strip("\n"))

		# For the last read of the last contig only: 
		nuc_covs = 0
		for pos in range(0, (len(contigs[contig_current])-1)):
			nuc_covs += sum(contigs[contig_current][pos].values())
		cov_average = float(nuc_covs)/float(contigs[contig_current][-1])
		contigs[contig_current].append(cov_average)

		real_insertions_coordinates = []
		real_insertions_motifs = []
		for ins in sorted(set(insertions[contig_current])):
			if contigs[contig_current][-2]-1 > ins[0]:
				cov_at_edges = float(sum(contigs[contig_current][ins[0]].values())+sum(contigs[contig_current][ins[0]+1].values()))/2
			else:
				cov_at_edges = float(sum(contigs[contig_current][ins[0]].values()))

			if cov_at_edges >= min_depth:
				if (insertions[contig_current].count(ins) >= (cov_at_edges-insertions[contig_current].count(ins))) and (insertions[contig_current].count(ins) >= cov_at_edges*cons_threshold):
					real_insertions_coordinates.append(ins[0])
					real_insertions_motifs.append(ins[1])
					print "Insertion detected in contig "+contig_current+", coverage at sides of insertion: "+str(cov_at_edges)+", coverage of the insertion: "+str(insertions[contig_current].count(ins))+", position/motif: "+str(ins)
		
		if real_insertions_coordinates != []:
			if sam_verify == True:
				print contig_current+" contains insertion(s), a separate SAM file will be additionally created for this contig."
				sam_file.append("@SQ"+"\t"+"SN:"+contig_current+"\t"+"LN:"+str(contigs[contig_current][-2]))
				for read in sam_reads:
					sam_file.append(read)
				outfile = open(outfolder+contig_current+"_to_"+prefix+".sam", "w")
				outfile.write("\n".join(sam_file)+"\n")
			insertions[contig_current] = [real_insertions_coordinates,real_insertions_motifs]
		else:
			del insertions[contig_current]

		print "Contig "+contig_current+" processed\n"

	mapfile.close()



	# Dictionary to translate IUPAC ambiguities
	amb = {"-A":"a",
		   "-C":"c",
		   "-G":"g",
		   "-N":"-",
		   "-T":"t",
		   "AC":"M",
		   "AG":"R",
		   "AN":"a",
		   "AT":"W",
		   "CG":"S",
		   "CN":"c",
		   "CT":"Y",
		   "GN":"g",
		   "GT":"K",
		   "NT":"t",
		   "-AC":"m",
		   "-AG":"r",
		   "-AN":"a",
		   "-AT":"w",
		   "-CG":"s",
		   "-CN":"c",
		   "-CT":"y",
		   "-GN":"g",
		   "-GT":"k",
		   "-NT":"t",
		   "ACT":"H",
		   "AGT":"D",
		   "CGT":"B",
		   "ACG":"V",
		   "-ACT":"h",
		   "-AGT":"d",
		   "-CGT":"b",
		   "-ACG":"v",       
		   "ACNT":"h",
		   "AGNT":"d",
		   "CGNT":"b",
		   "ACGN":"v",
		   "ACGT":"N"}



	# If a contig had zero reads mapped then set its coverage to 0.0 
	for contig in contigs:
		if type(contigs[contig][-2]) != int:
			contigs[contig].append(float(0))



	# Obtain sequence from the "contigs" dictionary
	fastas = {}
	for contig in contigs:
		for pos in range(0, contigs[contig][-2]):
			count_nucs = list(sorted(contigs[contig][pos].iteritems(), key=operator.itemgetter(1), reverse=True))
			cov_site = float(sum(contigs[contig][pos].values()))
			if cov_site < min_depth:
				if contig not in fastas:
					fastas[contig] = fill
				else:
					fastas[contig] += fill
			elif count_nucs[0][1] >= cons_threshold*cov_site:
				if contig not in fastas:
					fastas[contig] = count_nucs[0][0]
				else:
					fastas[contig] += count_nucs[0][0]
			elif count_nucs[0][1] + count_nucs[1][1] >= cons_threshold*cov_site:
				genotype = "".join(sorted([count_nucs[0][0],count_nucs[1][0]]))
				if contig not in fastas:
					fastas[contig] = amb[genotype]
				else:
					fastas[contig] += amb[genotype]
			elif count_nucs[0][1] + count_nucs[1][1] + count_nucs[2][1] >= cons_threshold*cov_site:
				genotype = "".join(sorted([count_nucs[0][0],count_nucs[1][0],count_nucs[2][0]]))
				if contig not in fastas:
					fastas[contig] = amb[genotype]
				else:
					fastas[contig] += amb[genotype]
			elif count_nucs[0][1] + count_nucs[1][1] + count_nucs[2][1] + count_nucs[3][1] >= cons_threshold*cov_site:
				genotype = "".join(sorted([count_nucs[0][0],count_nucs[1][0],count_nucs[2][0],count_nucs[3][0]]))
				if contig not in fastas:
					fastas[contig] = amb[genotype]
				else:
					fastas[contig] += amb[genotype]
			else:
				if contig not in fastas:
					fastas[contig] = "N"
				else:
					fastas[contig] += "N"



	# Add insertions for contigs with real insertions
	if insertions:
		for contig in insertions:
			if insertions[contig] != []:
				seq = fastas[contig]
				start = 0
				seqout = ""
				for i in range(0, len(insertions[contig][0])):
					seqout += seq[start:insertions[contig][0][i]]
					seqout += insertions[contig][1][i]
					start = insertions[contig][0][i]
				seqout += seq[insertions[contig][0][-1]:]
				fastas[contig] = seqout



	# Write fasta output files
	for contig in fastas:
		outfile = open(outfolder+prefix+"_to_"+contig+"_cons"+str(cons_threshold)+".fasta", "w")
		outfile.write(">"+prefix+" Mapped to: "+contig+", consensus threshold: "+str(cons_threshold)+", coverage: "+str(contigs[contig][-1])+"\n"+"\n".join(fastas[contig][i:i+nchar] for i in range(0, len(fastas[contig]), nchar))+"\n")
		outfile.close()



if __name__ == "__main__":
	main()
