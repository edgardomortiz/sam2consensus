#!/usr/bin/env python2
# -*- coding: utf-8 -*-


'''
+------------------------------------------------------------------+
| sam2consensus.py: extract the consensus sequence from a SAM file |
+------------------------------------------------------------------+

The program takes as input a SAM file (.sam or .sam.gz) resulting from mapping 
short reads to a reference (the reference sequences can correspond to separate 
genes for example), then it calculates the consensus sequence from the aligned
reads alone. A single or multiple consensus thresholds can be specified, the 
program also adds insertions, if many long insertions are expected, we recommend 
to perform indel ralignment before for optimal results. The consensus method is 
the one used by Geneious and described in detail in 
http://assets.geneious.com/manual/8.1/GeneiousManualse41.html

Regions with no coverage are filled with -s (or a different character if 
specified). Input SAM files don't need to be sorted. Original reference FASTAs 
are not necessary since the consensus is reference-free.

It will produce a FASTA file per reference containing as many sequences as 
thresholds were specified.
'''


__author__      = "Edgardo M. Ortiz"
__credits__     = "Deise J.P. Gonçalves"
__version__     = "2.1"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2019-01-11"



import sys
import os
import re
import operator
import argparse
import gzip
import math



def parsecigar(cigarstring, seq, pos_ref):
	'''
	Modifies a sequence according to its CIGAR string, modified for speed according to frequency of CIGAR tags

	:param cigarstring: cigar string, v3 or v4 are allowed
	:param seq: raw sequence
	:param pos_ref: position in the reference of the leftmost aligned nucleotide (0-based)
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
		if cigar[c]["type"] in ["=","X","M"]:
			seqout += seq[start:start + l]
			start += l
			start_ref += l
		elif cigar[c]["type"] in ["D","N","P"]:
			seqout += "-" * l
			start_ref += l
		elif cigar[c]["type"] == "I":
			insert.append((start_ref, seq[start:start+l]))
			start += l
		elif cigar[c]["type"] == "S":
			start += l
		elif cigar[c]["type"] == "H":
			continue
		else:
			print("SAM file probably contains unmapped reads")
	return seqout, insert



def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
		help="Name of the SAM file, SAM does not need to be sorted and can be compressed with gzip")
	parser.add_argument("-c", "--consensus-thresholds", action="store", dest="thresholds", type=str, default="0.25",
		help="List of consensus thresold(s) separated by commas, no spaces, example: -c 0.25,0.75,0.50, default=0.25")
	parser.add_argument("-n", action="store", dest="n", type=int, default=0,
		help="Split FASTA output sequences every n nucleotides, default=do not split sequence")
	parser.add_argument("-o", "--outfolder", action="store", dest="outfolder", default="./",
		help="Name of output folder, default=same folder as input")
	parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="",
		help="Prefix for output file name, default=input filename without .sam extension")
	parser.add_argument("-m", "--min-depth", action="store", dest="min_depth", type=int, default=1,
		help="Minimum read depth at each site to report the nucleotide in the consensus, default=1")
	parser.add_argument("-f", "--fill", action="store", dest="fill", default="-",
		help="Character for padding regions not covered in the reference, default= - (gap)")
	parser.add_argument("-d", "--maxdel", action="store", dest="maxdel", default=150,
		help="Ignore deletions longer than this value, default=150")
	args = parser.parse_args()



	filename = args.filename

	# Prepare the opener if the SAM file is compressed
	if filename.endswith(".gz"):
		opener = gzip.open
	else:
		opener = open

	# Parse list of consensus thresholds
	thresholds = args.thresholds.split(",")
	thresholds = [float(i) for i in thresholds]

	# Prefix will be input filename without extension
	if args.prefix == "":
		prefix = "".join(args.filename.split("/")[-1]).split(".")[0]
	else:
		prefix = args.prefix
	
	# Create output folder if it doesn't exist
	outfolder = args.outfolder.rstrip("/")
	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
	outfolder += "/"

	min_depth = args.min_depth

	nchar = args.n

	fill = args.fill

	maxdel = args.maxdel



	# Start processing SAM file
	print("\nProcessing file "+filename+":\n")

	# First process header of file to extract reference information

	header_length = 0

	with opener(filename) as mapfile:

		sequences  = {}             # Container of sequences per reference for consensus calculation
		coverages  = {}				# Container for coverage per site for each reference
		insertions = {}             # Container for insertions with coordinates per reference

		for line in mapfile:

			if line.startswith("@"):
				header_length += 1

				if line.startswith("@SQ"):

					# Obtain reference name and length, skipping description if present
					refname   = line.split("\t")[1].replace("SN:","").split()[0]
					reflength = int(line.split("\t")[2].replace("LN:",""))

					# Prepare empty dictionaries to parse aligned reads and populate them later
					sequences[refname]  = [{"-":0,"A":0,"C":0,"G":0,"N":0,"T":0} for nuc in range(reflength)]
					coverages[refname]  = [0]*reflength
					insertions[refname] = []

			else:
				break

	print("SAM header processed, "+str(len(sequences))+" references found.\n")
	mapfile.close()



	# Process the reads that were aligned to the references
	with opener(filename) as mapfile:
		# Set  read counter to show progress
		reads_total  = 0 - header_length
		reads_mapped = 0

		while 1:
			# Load large SAM file by chunks
			mapfile_chunk = mapfile.readlines(50000)
			if not mapfile_chunk:
				break

			for line in mapfile_chunk:

				# Skip header and unmapped reads
				reads_total += 1
				if line[0] != "@" and line.split("\t")[5] != "*":

					reads_mapped += 1
					sam_record  = line.split("\t")

					refname = sam_record[2].split()[0]        ## Reference name, skip description if mapper kept it
					pos_ref = int(sam_record[3])-1            ## Starting position in the reference, 0-based
					# cigar  = sam_record[5]                  ## CIGAR string fo the aligned read
					# seqraw = sam_record[9]                  ## Unaligned raw sequence of the read

					# Parse the CIGAR and obtain edited sequence and list of insertions
					seqout, insert = parsecigar(sam_record[5], sam_record[9], pos_ref)

					# Fill the nucleotides in their respective reference according to the alignment described by the
					# CIGAR string
					if seqout.count("-") <= maxdel:
						for nuc in seqout:
							sequences[refname][pos_ref][nuc] += 1
							pos_ref += 1
					else:
						for nuc in seqout:
							if nuc != "-":
								sequences[refname][pos_ref][nuc] += 1
							pos_ref += 1

					# Add insertions with coordinates to the dictionary of insertions per reference
					insertions[refname] += insert

				# Show progress
				if reads_total % 500000 == 0:
					print(str(reads_total)+" reads processed.")

	print("A total of "+str(reads_total)+" reads were processed, out of which, "+str(reads_mapped)+" reads were mapped.\n")
	mapfile.close()



	# Reformat dictionaries and get them ready for consensus calculation
	for refname in sequences:

		# Populate dictionary of coverage per site
		for pos in range(len(coverages[refname])):
			coverages[refname][pos] = sum(sequences[refname][pos].values())

			# Invert the dictionary at each position {"-":0,"A":0,"C":0,"G":0,"N":0,"T":0} to have counts as keys and a 
			# list of nucleotides with the same counts as values, skip nucs with count==0
			count_nucs = {}
			for key, value in sequences[refname][pos].iteritems():
				if value != 0:
					count_nucs.setdefault(value,[]).append(key)

			# Transform to a list of sorted tuples in reverse order
			count_nucs = list(count_nucs.iteritems())
			count_nucs.sort(reverse=True)

			# Multiply coverage by number of nucleotides in each count to restore the correct coverage
			count_nucs = [[i[0]*len(i[1]), i[1]] for i in count_nucs]
			sequences[refname][pos] = count_nucs



		# We will treat each insertion as a short alignment of motifs, for that we have to modify the raw counts of
		# motifs per coordinate and get them into the same shape as the dictionary that contains the nucleotide counts
		# for the rest of the sequence. Every column of each insertion will have a dictionary of the form:
		# {"-":0,"A":0,"C":0,"G":0,"N":0,"T":0} and will be parsed as the 'sequences' dictionary.

		# First step is to make a dictionary where the key is the inserted motif and the value is the count of the motif
		ins_tmp1 = {}
		if insertions[refname] != []:
			for ins in insertions[refname]:
				if ins[0] not in ins_tmp1:
					ins_tmp1[ins[0]] = {ins[1]:1}
				else:
					if ins[1] not in ins_tmp1[ins[0]]:
						ins_tmp1[ins[0]][ins[1]] = 1
					else:
						ins_tmp1[ins[0]][ins[1]] += 1

			# Now we create a new dictionary as long as the longest motif inserted for every coordinate where an
			# insertion was detected. Each position will contain a dictionary for counts as in the 'sequences'
			# dictionary: {"-":0,"A":0,"C":0,"G":0,"N":0,"T":0}
			ins_tmp2 = {}
			for pos in sorted(ins_tmp1):
				motif_max_len = 0
				for motif in ins_tmp1[pos]:
					motif_max_len = max(len(motif), motif_max_len)
				ins_tmp2[pos] = [{"-":0,"A":0,"C":0,"G":0,"N":0,"T":0} for i in range(motif_max_len)]
			
			# Now we track each nucleotide in each inserted motif and populate that new dictionary 'ins_tmp2'
			for pos in sorted(ins_tmp1):
				for motif in ins_tmp1[pos]:
					for col in range(len(motif)):
						ins_tmp2[pos][col][motif[col]] += ins_tmp1[pos][motif]

			# We are finally ready to reformat the nucleotide counts exactly as in the 'sequences' dictionary
			for pos in sorted(ins_tmp2):
				for col in range(len(ins_tmp2[pos])):

					# We need to add the count for "-" based on the coverage at the insertion position
					ins_tmp2[pos][col]["-"] = coverages[refname][pos] - sum(ins_tmp2[pos][col].values())

					# Invert the dictionary at each position {"A":0,"C":0,"T":0,"G":0,"-":0,"N":0} to have counts as 
					# keys and a list of nucleotides with the same counts as values, skip nucs with count==0
					count_nucs = {}
					for key, value in ins_tmp2[pos][col].iteritems():
						if value != 0:
							count_nucs.setdefault(value,[]).append(key)

					# Transform to a list of sorted tuples in reverse order
					count_nucs = list(count_nucs.iteritems())
					count_nucs.sort(reverse=True)

					# Multiply coverage by number of nucleotides in each count to restore the correct coverage
					count_nucs = [[i[0]*len(i[1]), i[1]] for i in count_nucs]
					ins_tmp2[pos][col] = count_nucs

			insertions[refname] = ins_tmp2



	# Dictionary to translate IUPAC ambiguities, lowercase letters are used when "-" or "N" were present for a position,
	# however, software like Genious for example are case insensitive and will imply ignore capitalization
	amb = {"-":"-", "A":"A", "C":"C", "G":"G", "N":"N", "T":"T",
		   "-A":"a", "-C":"c", "-G":"g", "-N":"n", "-T":"t",
		   "AC":"M", "AG":"R", "AN":"a", "AT":"W", "CG":"S",
		   "CN":"c", "CT":"Y", "GN":"g", "GT":"K", "NT":"t",
		   "-AC":"m", "-AG":"r", "-AN":"a", "-AT":"w", "-CG":"s",
		   "-CN":"c", "-CT":"y", "-GN":"g", "-GT":"k", "-NT":"t",
		   "ACG":"V", "ACN":"m", "ACT":"H", "AGN":"r", "AGT":"D",
		   "ANT":"w", "CGN":"s", "CGT":"B", "CNT":"y", "GNT":"k",
		   "-ACG":"v", "-ACN":"m", "-ACT":"h", "-AGN":"r", "-AGT":"d",
		   "-ANT":"w", "-CGN":"s", "-CGT":"b", "-CNT":"y", "-GNT":"k",
		   "ACGN":"v", "ACGT":"N", "ACNT":"h", "AGNT":"d", "CGNT":"b",
		   "-ACGN":"v", "-ACGT":"N", "-ACNT":"h", "-AGNT":"d", "-CGNT":"b",
		   "-ACGNT":"N"}



	# If a reference had zero reads mapped we need to erase it from the dictionaries
	sequences2erase = []
	for refname in coverages:
		if sum(coverages[refname]) == 0:
			sequences2erase.append(refname)
	for refname in sequences2erase:
		del sequences[refname]
		del insertions[refname]



	# Obtain consensus sequence(s) from the 'sequences' and 'insertions' dictionary
	fastas = {}
	for refname in sequences:

		for t in thresholds:
		
			fasta_seqout = ""
			fasta_header = ""
			sumcov = 0

			# Obtain sequence from the 'sequences' dictionary
			for pos in range(len(sequences[refname])):
				if sequences[refname][pos] != []:
					sumcov += coverages[refname][pos]
					if coverages[refname][pos] >= min_depth:
						nucs = []
						cov_nucs = 0
						for count in sequences[refname][pos]:
							if cov_nucs < t*coverages[refname][pos]:
								nucs     += count[1]
								cov_nucs += count[0]
							else:
								break
						fasta_seqout += amb["".join(sorted(nucs))]

						# Add insertions when applicable
						if refname in insertions:
							if pos in insertions[refname]:
								for col in range(len(insertions[refname][pos])):
									nucs = []
									cov_nucs = 0
									for count in insertions[refname][pos][col]:
										if cov_nucs < t*coverages[refname][pos]:
											nucs     += count[1]
											cov_nucs += count[0]
										else:
											break
									if amb["".join(sorted(nucs))] == "-":
										continue
									else:
										fasta_seqout += amb["".join(sorted(nucs))]
										sumcov += coverages[refname][pos]
					else:
						fasta_seqout += fill
				else:
					fasta_seqout += fill

			# Prepare headers for FASTA file
			# sequence name is: sammplename|consensus_threshold
			# sequence description is reference:refname coverage:XXX.XX length:XXXX consensus_threshold:XX%
			fasta_header = (">"+prefix+"|c"+str(int(t*100))+" reference:"+refname+
							" coverage:"+str(round(float(sumcov)/float(len(fasta_seqout)), 2))+
							" length:"+str(len(fasta_seqout.replace("-","")))+
							" consensus_threshold:"+str(int(t*100))+"%")

			# Add only if sequence is not empty
			if len(fasta_seqout.replace("-","")) > 0:
				if refname not in fastas:
					fastas[refname] = [[fasta_header,fasta_seqout]]
				else:
					fastas[refname].append([fasta_header,fasta_seqout])
			else:
				continue



	# Write output files, one FASTA file per reference, containing as many sequences as consensus thresholds were asked
	for reference in fastas:
		outnameprefix  = reference+"__"+prefix
		outfile = open(outfolder+outnameprefix+".fasta", "w")
		if nchar == 0:
			outfile.write("\n".join([i[0]+"\n"+i[1] for i in fastas[reference]])+"\n")
		else:
			outfile.write("\n".join([i[0]+"\n"+"\n".join([i[1][s:s+nchar] for s in range(0, len(i[1]), nchar)]) for i in fastas[reference]])+"\n")
		outfile.close()
		if len(thresholds) == 1:
			print("Consensus sequence at "+str(int(thresholds[0]*100))+"% saved for "+
				  reference+" in: "+outfolder+outnameprefix+".fasta")
		else:
			print("Consensus sequences at "+",".join([str(int(i*100))+"%"  for i in thresholds])+" saved for "+
				  reference+" in: "+outfolder+outnameprefix+".fasta")

	print("Done.\n")


if __name__ == "__main__":
	main()
