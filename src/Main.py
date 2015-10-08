'''
Created on Sep 23, 2011

@author: ksahlin
'''

import warnings
import os

def ReadInContigseqs(contigfile):
    cont_dict = {}
    k = 0
    temp = ''
    accession = ''
    for line in contigfile:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            cont_dict[accession] = ''
            k += 1
        elif line[0] == '>':
            cont_dict[accession] = temp
            temp = ''
            accession = line[1:].strip().split()[0]
        else:
            temp += line.strip()
    cont_dict[accession] = temp
    return(cont_dict)

def Main(contigfile_, tuple_of_bamfiles, tuple_of_means, tuple_of_thresholds, edge_support, read_len, cont_threshold,
          ratio, output_dest, std_dev, covcutoff, haplratio, haplthreshold, detect_haplotype, detect_duplicate, gff_file, fosmidpool, mapquality):

    from time import time
    import CreateGraph as CG
    import MakeScaffolds as MS
    import GenerateOutput as GO
    import Parameter
    from copy import deepcopy
    tot_start = time()
    F = [] #list of (ordered) lists of tuples containing (contig_name, direction, position, length, links). The tuple is a contig within a scaffold and the list of tuples is the scaffold.
    Scaffolds = {}     #scaffold dict with contig objects for easy fetching of all contigs in a scaffold
    n = len(tuple_of_bamfiles) # number of libraries we have
    param = Parameter.parameter() # object containing all parameters (user specified, defaulted and comuted along tha way.)
    param.scaffold_indexer = 1 # global indicator for scaffolds, used to index scaffolds when they are created
    param.map_quality = mapquality
    param.rel_weight = ratio
    Contigs = {} # contig dict that stores contig objects

    if not os.path.exists(output_dest):
        os.makedirs(output_dest)
    param.information_file = open(os.path.join(output_dest + 'Statistics.txt'), 'w')

    Information = param.information_file
    open(output_dest + '/haplotypes.fa', 'w')
    #Read in the sequences of the contigs in memory
    contigfile = open(contigfile_, 'r')
    C_dict = ReadInContigseqs(contigfile)
    #C_dict = {}
    param.gff_file = gff_file
#iterate over libraries
    param.first_lib = True
    for i in range(0, n):
        start = time()
        param.bamfile = tuple_of_bamfiles[i]
        param.mean_ins_size = tuple_of_means[i]
        param.ins_size_threshold = tuple_of_thresholds[i]
        param.edgesupport = edge_support[i]
        param.read_len = read_len[i]
        param.output_directory = output_dest
        param.std_dev_ins_size = std_dev[i]
        param.contig_threshold = cont_threshold[i]
        param.cov_cutoff = covcutoff[i]
        param.hapl_ratio = haplratio
        param.hapl_threshold = haplthreshold
        param.detect_haplotype = detect_haplotype
        param.detect_duplicate = detect_duplicate
        param.fosmidpool = fosmidpool
        print >> Information, '\nPASS ' + str(i + 1) + '\n\n'
        print 'Starting scaffolding with library: ', param.bamfile
        (G, Contigs, Scaffolds, F, param) = CG.PE(Contigs, Scaffolds, F, Information, output_dest, C_dict, param)      #Create graph, single out too short contigs/scaffolds and store them in F
        param.first_lib = False   #not the first lib any more
        if G == None:
            print '0 contigs/super-contigs passed the length criteria of this step. Exiting and printing results.. '
            break
        elapsed = time() - start
        print >> Information, 'Time elapsed for creating graph, iteration ' + str(i) + ': ' + str(elapsed) + '\n'
        start = time()
        (Contigs, Scaffolds, F, param) = MS.Algorithm(G, Contigs, Scaffolds, F, Information, C_dict, param)   # Make scaffolds, store the complex areas (consisting of contig/scaffold) in F, store the created scaffolds in Scaffolds, update Contigs
        elapsed = time() - start
        print >> Information, 'Time elapsed for making scaffolds, iteration ' + str(i) + ': ' + str(elapsed) + '\n'

        print 'Writing out scaffolding results for step', i + 1, ' ...'
        Scaffolds_copy = deepcopy(Scaffolds)
        Contigs_copy = deepcopy(Contigs)
        F_copy = deepcopy(F)
        for scaffold_ in Scaffolds_copy.keys(): #iterate over keys in hash, so that we can remove keys while iterating over it
            ###  Go to function and print to F
            ### Remove Scaf_obj from Scaffolds and Contig_obj from contigs
            S_obj = Scaffolds_copy[scaffold_]
            list_of_contigs = S_obj.contigs   #list of contig objects contained in scaffold object
            Contigs_copy, F_copy = GO.WriteToF(F_copy, Contigs_copy, list_of_contigs)
            del Scaffolds_copy[scaffold_]
        #print F
        GO.PrintOutput(F_copy, C_dict, Information, output_dest, param, i + 1)

    ### Calculate stats for last scaffolding step    
    scaf_lengths = [Scaffolds[scaffold_].s_length for scaffold_ in Scaffolds.keys()]
    sorted_lengths = sorted(scaf_lengths, reverse=True)
    N50, L50 = CG.CalculateStats(sorted_lengths, param)
    param.current_L50 = L50
    param.current_N50 = N50
#    ### Call a print scaffolds function here for remaining scaffolds that has "passed" all library levels
#    for scaffold_ in Scaffolds.keys(): #iterate over keys in hash, so that we can remove keys while iterating over it
#        ###  Go to function and print to F
#        ### Remove Scaf_obj from Scaffolds and Contig_obj from contigs
#        S_obj=Scaffolds[scaffold_]
#        list_of_contigs=S_obj.contigs   #list of contig objects contained in scaffold object
#        Contigs, F = GO.WriteToF(F,Contigs,list_of_contigs)
#        del Scaffolds[scaffold_]
#    #print F
#    GO.PrintOutput(F,C_dict,Information,output_dest)

    elapsed = time() - tot_start
    print >> Information, 'Total time for scaffolding: ' + str(elapsed) + '\n'
    print 'Finished\n\n '





if __name__ == '__main__':
    import sys
    nr_of_libs = int(sys.argv[1])
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-s", "--stddev", dest="stddev", default=(None,) * nr_of_libs,
                      help="estimated standard deviation of libraries", nargs=nr_of_libs, type="int")

    parser.add_option("-o", "--output", dest="output", default='../',
                      help="The output location", type="string")

    parser.add_option("-f", "--bamfile", dest="bamfiles",
                      help="Name of bamfile", nargs=nr_of_libs, type="string")

    parser.add_option("-r", dest="readlen", nargs=nr_of_libs, default=(None,) * nr_of_libs,
                      help="Mean read length ", type="int")

    parser.add_option("-m", dest="mean", default=(None,) * nr_of_libs,
                      help="mean of insert library", nargs=nr_of_libs, type="int")

    parser.add_option("-w", dest="relweight", default=3,
                      help="treshold value for the relative weight of an edge", type="int")

    parser.add_option("-T", dest="threshold", nargs=nr_of_libs, default=(20000,) * nr_of_libs,
                      help="treshold value ", type="int")

    parser.add_option("-c", dest="contigfile",
                      help="file of contigs", type="string")

    parser.add_option("-e", dest="edgesupport", nargs=nr_of_libs, default=(5,) * nr_of_libs,
                      help="treshold value for the least nr of links that is needed to create an edge. Default for all libs: 5 ", type="int")
    parser.add_option("-k", "--minsize", dest="minsize", type="int", nargs=nr_of_libs, default=(500,) * nr_of_libs,
                  help="contig size threshold for the library (all contigs below this size is discarded from the scaffolding). Default: Set to same as -T parameter")

    parser.add_option("-z", "--covcutoff", dest="covcutoff", default=(None,) * nr_of_libs,
                      help="User specified coverage cutoff", nargs=nr_of_libs, type="int")

    parser.add_option("-a", "--haplratio", dest="haplratio", default=1.3,
                      help="Maximum length ratio for merging of haplotypic regions", type="float")

    parser.add_option("-b", "--haplthreshold", dest="haplthreshold", default=3,
                      help="Nr of standard deviations over mean/2 of coverage to allow for clasification of haplotype", type="int")

    parser.add_option("-g", "--haplotype", dest="haplotype", default=0,
                      help="Haplotype detection function, default = on", type="int")

    parser.add_option("-d", "--duplicate", dest="duplicate", default=1,
                      help="Sequencing duplicates detection, default = off", type="int")

    parser.add_option("--mapq", "--lowestmapquality", dest="mapquality", default=10,
                      help="Lowest mapping quality allowed in order to use the read. This value is compared to the mapping quality column in the BAM file. ", type="int")

#TEMPORARY parameters here, remove after spruce assembly
    parser.add_option("-t", dest="transcriptalignfile",
                      help="file of contigs", type="string")

    parser.add_option("-p", dest="fosmidpool", default=None,
                      help="""Specify that data comes from a fosmid pool. This parameter is sets the number of
                      links that the second strongest link must have in order to not do any scaffolding
                      with the region.""", type="int")
    (options, args) = parser.parse_args()

    #print options.transcriptalignfile
    if nr_of_libs == 1:
        options.transcriptalignfile = (options.transcriptalignfile,)
        options.bamfiles = (options.bamfiles,)
        #options.threshold=(options.threshold,)
        options.covcutoff = (options.covcutoff,)
        #options.mean=(options.mean,)
        #options.stddev=(options.stddev,)
        if type(options.stddev) == int:
            options.stddev = (options.stddev,)
        if type(options.mean) == int:
            options.mean = (options.mean,)
        if type(options.threshold) == int:
            options.threshold = (options.threshold,)
        if type(options.readlen) == int:
            options.readlen = (options.readlen,)
        if type(options.minsize) == int:
            options.minsize = (options.minsize,)
        if type(options.edgesupport) == int:
            options.edgesupport = (options.edgesupport,)

    if (options.mean[0] and not options.stddev[0]) or (not options.mean[0] and options.stddev[0]):
        parser.error("either both or none of -m and -s is required")
    if (options.threshold[0] and not options.minsize[0]) or (not options.threshold[0] and options.minsize[0]):
        parser.error("either both or none of -T and -t is required")
    if not options.contigfile:
        parser.error("parameter -c (a fasta contig file) is required")
    if not options.bamfiles:
        parser.error("parameter -f (BAM files) is required")
    if not options.haplotype and (options.haplratio or options.haplthreshold):
        warnings.warn('parameter -g (treating haplotypic regions) inactivated, parameters -a and -b will not have any effect if specified. ')

    #check that bam files exists
    for file_ in options.bamfiles:
        try:
            open(file_)
        except IOError as e:
            sys.exit("couldn't find BAM file: " + file_ + " check that the path is correct and that the file exists")
        try:
            open(file_ + '.bai')
        except IOError as e:
            print "couldn't find index file: ", file_ + '.bai', " check that the path is correct and that the file exists"
            sys.exit(0)
    #check that contig files exists
    try:
        open(options.contigfile)
    except IOError as e:
        sys.exit("couldn't open contig file " + options.contigfile + " check that the path is correct and that the file exists")


    Main(options.contigfile, options.bamfiles, options.mean, options.threshold, options.edgesupport, options.readlen,
          options.minsize, options.relweight, options.output, options.stddev, options.covcutoff, options.haplratio,
          options.haplthreshold, options.haplotype, options.duplicate, options.transcriptalignfile, options.fosmidpool,
          options.mapquality)


