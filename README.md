BESST_RNA
=========

Scaffolding of genomic assemblies with RNA seq data.


INPUT:
-------
Required arguments:

* -c <str> Path to a contig file. 

*  -f <str> Path to bamfile of RNA-seq data mapped to assembly. 

* -o <str> Path to location for the output.

Optional arguments

* -e <int> The least amount of witness links that is needed to create a link edge in graph (one for each library). 

* -w <int> real weight parameter. Given multiple edges from a given contig(-end) inference on which one to be chosen needs to be done. -w is the minimum threshold of the relative weight for the dominating edge (most amount of links) to the second most dominating edge for a given contig-end. If weight is less than specified -w, no edge will be chosen to extend this contig. Default value is 3.

* -T <int> Upper threshold on the distance between RNA reads. RNA reads with distance over this threshold is not considered to be mapped correctly. e.g. could be set to the upper estimated intron size.

* -k <int> Minimum size of contig to be included in the scaffolding.

* -d <bool> Check for sequencing duplicates and count only one of them (when computing nr of links) if they are occurring. (default on = 1).

* -z <int> Coverage cutoff for repeat classification > Should be set high (i.e repeat detection turned of) for RNA scaffolding since repeat classification is not trusted with this type of data. [ e.g. -z 1000 says that contigs with coverage over 100 will be discarded from scaffolding.]

* --mapq <int> Lowest mapping quality allowed in order to use the read. This value is compared to the mapping quality column in the BAM file. Default = 10.

EXAMPLE RUN:

From /<path to>/BESST_RNA/src/ , run

$ python Main.py 1 -c /path/to/contigfile.fa -f /path/to/bamfile -o /path/to/output	-e 3  -T 20000 -k 500  -d 1  -z 1000 



NOTE:
-------

Mapping reads: BESST assumes "PE orientation" the following mapping of RNA reads:

PE: 
   s                    t
   ------>      <-------









