# sam2consensus
Get the consensus sequences for reads mapped to a reference made of multiple separate genes.

## _Brief description_
The program takes as input a SAM file resulting from mapping short reads to a collection of
gene sequences as reference, then it calculates the consensus sequence per gene without
considering the reference. It adds insertions and can take a custom consensus threshold,
the consensus method is the same as the one described for [Geneious](http://assets.geneious.com/manual/8.1/GeneiousManualse41.html). Regions that were not covered will be filled with Ns.

Input SAM files have to be sorted an can be compressed with `gzip` (BAM is not supported). The processing is faster if the SAM file contains only mapped reads. Original reference FASTAs are not necessary.

It will produce a FASTA sequence per gene, and in case the gene has insertions it will also create
a separate SAM file just for the particular gene for verification purposes.

## _Usage_
Just type `python sam2consensus.py -h` to show the help of the program:
```
usage: sam2consensus.py [-h] -i FILENAME [-c CONS_THRESHOLD] [-m MIN_DEPTH]
                        [-o OUTFOLDER] [-p PREFIX] [-f FILL] [-s]

Calculates the consensus sequence from reads aligned to a multi-gene fasta
reference

optional arguments:
  -h, --help            show this help message and exit
  -i FILENAME, --input FILENAME
                        Name of the SAM file, SAM must be sorted before!
  -c CONS_THRESHOLD, --consensus-threshold CONS_THRESHOLD
                        Consensus threshold sensu Geneious, default=0.25
  -m MIN_DEPTH, --min-depth MIN_DEPTH
                        Minimum read depth at each site to report the
                        nucleotide in the consensus, default=5
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Name of output folder, default=same folder as input
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name, default=input filename
                        without .sam extension
  -f FILL, --fill FILL  Character for padding regions not covered in the
                        reference, default= - (gap)
  -s, --sam-verify      Enable creation of individual SAM files per gene for
                        verification of insertion, disabled by default
```

## _Examples_
_Example 1:_ Using the default consensus threshold of 0.25:
```bash
python sam2consensus.py -i myfile.sam
```

_Example 2:_ Using a custom consensus threshold of 0.50, and increasing the minimum depth of coverage to 10:
```bash
python sam2consensus.py -i myfile.sam -c 0.5 -m 10
```

_Example 3:_ Using a custom consensus threshold of 0.75 and specifying a different folder for output:
```bash
python sam2consensus.py -i myfile.sam -c 0.75 -o ./outfiles/
```

_Example 4:_ Specifying a prefix "sample1" and filling gaps with Ns:
```bash
python sam2consensus.py -i myfile.sam -pre sample1 -f N
```


## _Our pipeline for obtaining the SAM file_

## _Credits_
- Code: [Edgardo M. Ortiz](mailto:e.ortiz.v@gmail.com)
- Data and testing: [Deise J. P. Gonçalves](mailto:deisejpg@gmail.com)
