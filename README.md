# sam2consensus
Get the consensus sequences for short sequencing reads mapped to a reference. Version 2.0

## _Brief description_
The program takes as input a SAM file (.sam or .sam.gz) resulting from mapping short reads to a reference (the reference
sequences can correspond to separate genes for example), then it calculates the consensus sequence from the aligned
reads alone. If you have BAM files you will need to convert them first with `samtools`. A single or multiple consensus thresholds can be specified, the program also adds insertions, if many long
long insertions are expected, we recommend to perform indel ralignment before for optimal results. The consensus method
is the one used by Geneious and described in detail in http://assets.geneious.com/manual/8.1/GeneiousManualse41.html

Regions with no coverage are filled with -s (or a different character if specified). Input SAM files don't need to be
sorted. Original reference FASTAs are not necessary since the consensus is reference-free.

It will produce a FASTA file per reference containing as many sequences as thresholds were specified.

## _Usage_
Just type `python sam2consensus.py -h` to show the help of the program:
```
usage: sam2consensus.py [-h] -i FILENAME [-c THRESHOLDS] [-n N] [-o OUTFOLDER]
                        [-p PREFIX] [-m MIN_DEPTH] [-f FILL] [-d MAXDEL]

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

optional arguments:
  -h, --help            show this help message and exit
  -i FILENAME, --input FILENAME
                        Name of the SAM file, SAM does not need to be sorted
                        and can be compressed with gzip
  -c THRESHOLDS, --consensus-thresholds THRESHOLDS
                        List of consensus thresold(s) separated by commas, no
                        spaces, example: -c 0.25,0.75,0.50, default=0.25
  -n N                  Split FASTA output sequences every n nucleotides,
                        default=do not split sequence
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Name of output folder, default=same folder as input
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name, default=input filename
                        without .sam extension
  -m MIN_DEPTH, --min-depth MIN_DEPTH
                        Minimum read depth at each site to report the
                        nucleotide in the consensus, default=1
  -f FILL, --fill FILL  Character for padding regions not covered in the
                        reference, default= - (gap)
  -d MAXDEL, --maxdel MAXDEL
                        Ignore deletions longer than this value, default=150
```

## _Examples_
In the following examples, if you added the script to the `PATH` you can omit `python` from the commands.

_Example 1:_ Using the default consensus threshold of 0.25:
```bash
python sam2consensus.py -i myfile.sam
```

_Example 2:_ Using multiple consensus thresholds of 0.25, 0.50, and 0.75 and specifying an output folder:
```bash
python sam2consensus.py -i myfile.sam -c 0.25,0.50,0.75 -o results
```

_Example 3:_ Using a custom consensus threshold of 0.75 and changing padding character for uncovered regions to "N":
```bash
python sam2consensus.py -i myfile.sam -c 0.75 -f N
```

_Example 4:_ Specifying a prefix "sample1" and increasing minimum mapping depth to 10:
```bash
python sam2consensus.py -i myfile.sam -pre sample1 -m 10
```

## _Credits_
- Code: [Edgardo M. Ortiz](mailto:e.ortiz.v@gmail.com)
- Data and testing: [Deise J. P. Gonçalves](mailto:deisejpg@gmail.com)
