# sam2consensus
Get the consensus sequences for reads mapped to a reference made of multiple separate genes.

## Brief description
The program takes as input a SAM file resulting from mapping short reads to a collection of
gene sequences as reference, then it calculates the consensus sequence per gene without
considering the reference. It adds insertions and can take a custom consensus threshold,
the consensus method the same as in Geneious (http://assets.geneious.com/manual/8.1/GeneiousManualse41.html)

Input SAM files have to be sorted, contain only mapped reads (preferably), and have
CIGAR strings in SAM v.1.3. Original reference FASTAs are not necessary.

It will produce a FASTA sequence per gene, and in case the gene has insertions it will also create
a separate SAM file just for the particular gene for verification purposes.

## Usage
Simply specify the name of the SAM input file and optionally the consensus thresold in decimals:

_Example 1:_ Using the default consensus threshold of 0.25:
python sam2consensus.py myfile.sam

_Example 2:_ Using a custom consensus threshold of 0.50:
python sam2consensus.py myfile.sam 0.5

## Our pipeline for obtaining the SAM file
