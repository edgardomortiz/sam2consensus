# sam2consensus
Get the consensus sequences for reads mapped to a reference made of multiple separate genes.

## _Brief description_
The program takes as input a SAM file resulting from mapping short reads to a collection of
gene sequences as reference, then it calculates the consensus sequence per gene without
considering the reference. It adds insertions and can take a custom consensus threshold,
the consensus method is the same as the one described for [Geneious](http://assets.geneious.com/manual/8.1/GeneiousManualse41.html). Regions that were not covered will be filled with Ns.

Input SAM files have to be sorted and contain only mapped reads (preferably). Original reference FASTAs are not necessary.

It will produce a FASTA sequence per gene, and in case the gene has insertions it will also create
a separate SAM file just for the particular gene for verification purposes.

## _Usage_
Simply specify the name of the SAM input file and optionally the consensus threshold in decimals:

_Example 1:_ Using the default consensus threshold of 0.25:
```bash
python sam2consensus.py myfile.sam
```

_Example 2:_ Using a custom consensus threshold of 0.50:
```bash
python sam2consensus.py myfile.sam 0.5
```

_Example 3:_ Using a custom consensus threshold of 0.75 and specifying a different folder for output:
```bash
python sam2consensus.py myfile.sam 0.75 ./outfiles/
```

## _Our pipeline for obtaining the SAM file_
