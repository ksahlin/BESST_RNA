BESST_RNA
=========

Scaffolding of genomic assemblies with RNA seq data


BESST
======
 
INPUT:
-------
Required arguments:

* -c <path to a contig file>  

*  -f < path to bamfiles>  (increasing order of insert size)

* -o <path to location for the output>

* -e <The least amount of witness links that is needed to create a link edge in graph (one for each library)> (integer number) 

* -T <Upper threshold on the distance between RNA reads. RNA reads with distance over this threshold is not considered to be mapped correctly. e.g. could be set to the upper estimated intron size.

* -k  <Minimum size of contig to be included in the scaffolding> 

* -d <check for sequencing duplicates and count only one of them (when computing nr of links) if they are occurring> (default on = 1). <0 or 1> 

* -z <Coverage cutoff for repeat classification> Should be set high (i.e repeat detection turned of) for RNA scaffolding since repeat classification is not trusted with this type of data. [ e.g. -z 100 says that contigs with coverage over 100 will be discarded from scaffolding.]



EXAMPLE RUN:

runBESST -c /path/to/contigfile.fa -f /path/to/file1.bam /path/to/file2.bam -o /path/to/output	-e 3  -T 20000 -k 500  -d 1  -z 1000        -o /path/to/output

Optional arguments:

The following arguments are computed internally by BESST. It is however good to specify mean and standard deviation if your assembly is very fragmented compared to the library insert size (not enough large contains to compute library statistics on).


NOTE:
-------

Mapping reads: BESST assumes "PE orientation" the following mapping of RNA reads:

PE: 
   s                    t
   ------>      <-------









