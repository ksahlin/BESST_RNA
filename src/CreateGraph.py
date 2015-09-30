'''
Created on Sep 23, 2011

@author: ksahlin
'''

import pysam
import Contig, Scaffold, Parameter
import GenerateOutput as GO
import networkx as nx
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    pass

from time import time

def PE(Contigs, Scaffolds, F, Information, output_dest, C_dict, param):
    G = nx.Graph()
    print 'Parsing BAM file...'
    #informative_pair={81:(False,True),97:(True,False),113:(False,False),65:(True,True)}
    #I switched to look at mates instead since BWA can give false flag combinations for
    # read-mate when read is mapped but not mate eg 97-149 81-165. But the reverse
    #does not happen.
    #informative_pair={161:(True,False),145:(False,True),129:(True,True),177:(False,False)} #,131:(True,True),179:(False,False)} #147:(False,True),163:(True,False),
    with pysam.Samfile(param.bamfile, 'rb') as bam_file:    #once real data, change to 'rb', simulated files are on SAM format

        #Get parameters -r, -m, -s, -T, -t for library
        print 'Computing parameters not set by user...'
        GetParams(bam_file, param, Scaffolds, C_dict, F, Contigs)

        #Clean contig_library
        singeled_out = 0
        if param.first_lib:
            cont_lengths = bam_file.lengths
            cont_lengths = [int(nr) for nr in cont_lengths]  #convert long to int object
            cont_names = bam_file.references

            #Calculate NG50 and LG 50
            param.tot_assembly_length = sum(cont_lengths)
            sorted_lengths = sorted(cont_lengths, reverse=True)
            NG50, LG50 = CalculateStats(sorted_lengths, param)
            param.current_LG50 = LG50
            param.current_NG50 = NG50


####### WHEN ADDING SHORTER CONTIGS NOT INCLUDED IN THE SCAFFOLDING, 
####### WE NEED TO ALSO INITIALIZE OBJECTS FOR THESE, THIS SHOULD BE DONE SOMEWHERE HERE
            for i in range(0, len(cont_names)):
                if cont_lengths[i] >= param.contig_threshold:
                    C = Contig.contig(cont_names[i])   # Create object contig
                    C.length = cont_lengths[i]
                    scaf_length = C.length        # Initially, scaffold consists of only this contig    
                    C.direction = True              # always in same direction first, False=reverse
                    C.position = 0                  #position always 0
                    C.links = {}
                    Contigs[C.name] = C              # Create a dict with name as key and the object container as value
                    S = Scaffold.scaffold(param.scaffold_indexer, [C], scaf_length, {}, {})  # Create object scaffold
                    Scaffolds[S.name] = S
                    C.scaffold = S.name
                    param.scaffold_indexer += 1
                else:
                    singeled_out += 1
                    F.append([(cont_names[i], True, 0, cont_lengths[i], {})])   #list of (contig_name, pos_direction, position,length)
            print >> Information, 'Nr of contigs that was singeled out due to length constraints ' + str(singeled_out)
        else:
                #Clean contig_library/scaffold_library
            scaf_lengths = [Scaffolds[scaffold_].s_length for scaffold_ in Scaffolds.keys()]
            sorted_lengths = sorted(scaf_lengths, reverse=True)
            NG50, LG50 = CalculateStats(sorted_lengths, param)
            param.current_LG50 = LG50
            param.current_NG50 = NG50
            for scaffold_ in Scaffolds.keys(): #iterate over keys in hash, so that we can remove keys while iterating over it
                if Scaffolds[scaffold_].s_length < param.contig_threshold:
                    ###  Go to function and print to F
                    ### Remove Scaf_obj from Scaffolds and Contig_obj from contigs
                    S_obj = Scaffolds[scaffold_]
                    list_of_contigs = S_obj.contigs   #list of contig objects contained in scaffold object
                    Contigs, F = GO.WriteToF(F, Contigs, list_of_contigs)  #Don't worry, the contig objects are removed in WriteTOF function
                    del Scaffolds[scaffold_]
                    singeled_out += 1
            print >> Information, 'Nr of contigs/scaffolds that was singeled out due to length constraints ' + str(singeled_out)


        #Create "node graph" of contigs (that passed the length criteria). Each having a left and right node
        print 'Nr of contigs/scaffolds included in scaffolding: ' + str(len(Scaffolds))#,Scaffolds.keys()
        if len(Scaffolds) == 0:
            return(None, Contigs, Scaffolds, F, param)
        cnt = 0
        tot_start = time()
        start1 = time()
        for scaffold_ in Scaffolds:
            G.add_edge((scaffold_, 'L'), (scaffold_, 'R'), nr_links=None)    #this is a scaffold object but can be both a single contig or a scaffold.
            Scaffolds[ scaffold_ ].scaffold_left_nbrs = {}
            Scaffolds[ scaffold_ ].scaffold_right_nbrs = {}
            if cnt % 100000 == 0 and cnt > 0:
                elapsed = time() - start1
                print >> Information, 'Total nr of keys added: ', cnt, 'Time for adding last 100 000 keys: ', elapsed
                start1 = time()
            cnt += 1
        print 'Total time elapsed: ', time() - tot_start
        # Create the link edges in the graph by fetching info from bam file

        cont_aligned_len = {}
        for contig in Contigs:
            cont_aligned_len[contig] = [0, Contigs[contig].length]

        count = 0
        non_unique = 0
        non_unique_for_scaf = 0
        nr_of_duplicates = 0
        prev_obs1 = -1
        prev_obs2 = -1
        reads_with_too_long_insert = 0
        #fishy_reads = {}
        for alignedread in bam_file:
            try: #check that read is aligned OBS: not with is_unmapped since this flag is fishy for e.g. BWA
                contig1 = bam_file.getrname(alignedread.rname)
                contig2 = bam_file.getrname(alignedread.mrnm)
            except ValueError:
                continue
            #contig1=bam_file.getrname(alignedread.rname)
            ## add to coverage computation if contig is still in the list of considered contigs
            try:
                cont_aligned_len[contig1][0] += alignedread.rlen
            except KeyError:
                pass
########## CREATE EDGES IN SCAFFOLD GRAPH ##########

            if contig1 != contig2 and alignedread.is_read2:
                #check how many non unique reads out of the useful ones (mapping to two different contigs)
                #This only works for BWA!! implement for other aligners as well
                if alignedread.mapq == 0:
                    non_unique += 1
                    #print contig1,contig2
                if contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold and alignedread.mapq > param.map_quality: # and alignedread.tags[0][1] == 'U':
                    #if alignedread.tags[0][1] != 'U':
                    #    non_unique_for_scaf += 1
                    if alignedread.mapq == 0:
                        non_unique_for_scaf += 1
                    count += 1
                    #(read_dir,mate_dir)=informative_pair[flag_type]
                    (read_dir, mate_dir) = (not alignedread.is_reverse, not alignedread.mate_is_reverse)
                    scaf1 = Contigs[contig1].scaffold
                    scaf2 = Contigs[contig2].scaffold
                    #Calculate actual position on scaffold here
                    #position1 cont/scaf1
                    cont_dir1 = Contigs[contig1].direction  #if pos : L if neg: R
                    cont1_pos = Contigs[contig1].position
                    readpos = alignedread.pos
                    cont1_len = Contigs[contig1].length
                    s1len = Scaffolds[scaf1].s_length
                    #position1 cont1/scaf1                        
                    cont_dir2 = Contigs[contig2].direction
                    cont2_pos = Contigs[contig2].position
                    matepos = alignedread.mpos
                    cont2_len = Contigs[contig2].length
                    s2len = Scaffolds[scaf2].s_length
                    (obs1, obs2, scaf_side1, scaf_side2) = PosDirCalculatorPE(cont_dir1, read_dir, cont1_pos, readpos, s1len, cont1_len, cont_dir2, mate_dir, cont2_pos, matepos, s2len, cont2_len, param.read_len)
                    if obs1 == prev_obs1 and obs2 == prev_obs2:
                        nr_of_duplicates += 1
                        if param.detect_duplicate:
                            continue

                    if obs1 + obs2 < param.ins_size_threshold:
#                        if obs1 == 3 or obs2 ==3:
#                            print alignedread.pos,alignedread.mpos, contig1, contig2, scaf1, scaf2, s1len,s2len
                        if scaf_side1 == 'R':
                            if (scaf2, scaf_side2) in Scaffolds[scaf1].right_nbrs_obs:
                                if obs1 < Scaffolds[scaf1].right_nbrs_obs[(scaf2, scaf_side2)]:
                                    Scaffolds[scaf1].right_nbrs_obs[(scaf2, scaf_side2)] = obs1
                            else:
                                Scaffolds[scaf1].right_nbrs_obs[(scaf2, scaf_side2)] = obs1
                        if scaf_side1 == 'L':
                            if (scaf2, scaf_side2) in Scaffolds[scaf1].left_nbrs_obs:
                                if obs1 < Scaffolds[scaf1].left_nbrs_obs[(scaf2, scaf_side2)]:
                                    Scaffolds[scaf1].left_nbrs_obs[(scaf2, scaf_side2)] = obs1
                            else:
                                Scaffolds[scaf1].left_nbrs_obs[(scaf2, scaf_side2)] = obs1
                        if scaf_side2 == 'R':
                            if (scaf1, scaf_side1) in Scaffolds[scaf2].right_nbrs_obs:
                                if obs2 < Scaffolds[scaf2].right_nbrs_obs[(scaf1, scaf_side1)]:
                                    Scaffolds[scaf2].right_nbrs_obs[(scaf1, scaf_side1)] = obs2
                            else:
                                Scaffolds[scaf2].right_nbrs_obs[(scaf1, scaf_side1)] = obs2
                        if scaf_side2 == 'L':
                            if (scaf1, scaf_side1) in Scaffolds[scaf2].left_nbrs_obs:
                                if obs2 < Scaffolds[scaf2].left_nbrs_obs[(scaf1, scaf_side1)]:
                                    Scaffolds[scaf2].left_nbrs_obs[(scaf1, scaf_side1)] = obs2
                            else:
                                Scaffolds[scaf2].left_nbrs_obs[(scaf1, scaf_side1)] = obs2

                        if (scaf2, scaf_side2) not in G[(scaf1, scaf_side1)]:
                            G.add_edge((scaf2, scaf_side2), (scaf1, scaf_side1), nr_links=1, gap_dist=obs1 + obs2)
                        else:
                            G.edge[(scaf1, scaf_side1)][(scaf2, scaf_side2)]['nr_links'] += 1
                            G.edge[(scaf1, scaf_side1)][(scaf2, scaf_side2)]['gap_dist'] += obs1 + obs2
                    else:
                        reads_with_too_long_insert += 1
                        #fishy_reads[alignedread.qname[:-1]]=[contig2,alignedread.is_read2]
                        ## add to haplotype graph here!!

                    prev_obs1 = obs1
                    prev_obs2 = obs2

                elif contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold:
########################Use to validate scaffold in previous step here ############
                    pass
#        print 'NR OF FISHY EDGES: ', len(fishy_reads)
        print 'USEFUL READS (reads mapping to different contigs): ', count
    #print 'Non unique portion out of "USEFUL READS"  (filtered out from scaffolding): ', non_unique
        #print 'Non unique used for scaf: ', non_unique_for_scaf
        print 'Reads with too large insert size from "USEFUL READS" (filtered out): ', reads_with_too_long_insert
        if param.detect_duplicate:
            print 'Number of duplicated reads indicated and removed: ', nr_of_duplicates

    ##### Calc coverage for all contigs with current lib here #####
        sum_x = 0
        sum_x_sq = 0
        n = 0
        for contig in cont_aligned_len:
            cont_coverage = cont_aligned_len[contig][0] / float(cont_aligned_len[contig][1])
                #print key, cont_aligned_len[key]/float(cont_lengths[i])
            try:
                Contigs[contig].coverage = cont_coverage
            except KeyError:
                pass
            sum_x += cont_coverage
            sum_x_sq += cont_coverage ** 2
            n += 1

        mean_cov, std_dev_cov = CalculateMeanCoverage(Contigs, param.first_lib, output_dest, param.bamfile)
        param.mean_coverage = mean_cov
        param.std_dev_coverage = std_dev_cov


    return(G, Contigs, Scaffolds, F, param)


def PreFilterEdges(G, Scaffolds, param):
    #### Pre filtering of edges here #### 
    pre_filtered = 0
    for node in G.nodes():
        min_cov = 1000000
        for contig in Scaffolds[ node[0] ].contigs:
            if contig.coverage and contig.coverage < min_cov:
                min_cov = contig.coverage
        k = 2 * param.read_len / float(min_cov)
        if node[1] == 'R':
            for nbr in G.neighbors(node): #for nbr in Scaffolds[node[0] ].right_nbrs_obs:
                #lowest k from node or nbr determines k
                if nbr[0] != node[0]:
                    try:
                        for contig in Scaffolds[ nbr[0] ].contigs:
                            if contig.coverage and contig.coverage < min_cov:
                                min_cov = contig.coverage
                        k = 2 * param.read_len / float(min_cov)
                    except KeyError:
                        continue
                    lower_bound = Scaffolds[node[0] ].right_nbrs_obs[nbr]
                    #print lower_bound
                    if lower_bound > param.read_len + 20 * k:
#                        for cont_obj in Scaffolds[node[0] ].contigs:
#                            print cont_obj.name
#                        try:
#                            nr_links = G[node][nbr]['nr_links'] 
#                            print 'links: ',nr_links
#                        except KeyError:
#                            try:
#                                nr_links = G[node][nbr]['nr_links']
#                                print 'links: ', nr_links
#                            except KeyError:
#                                print 'a repeat removed'
                        G.remove_edge(node, nbr)
                        pre_filtered += 1

        else:
            for nbr in G.neighbors(node): #Scaffolds[node[0] ].left_nbrs_obs:
                #lowest k from node or nbr determines k
                if nbr[0] != node[0]:

                    try:
                        for contig in Scaffolds[ nbr[0] ].contigs:
                            if contig.coverage and contig.coverage < min_cov:
                                min_cov = contig.coverage
                        k = 2 * param.read_len / float(min_cov)
                    except KeyError:
                        continue
                    lower_bound = Scaffolds[node[0] ].left_nbrs_obs[nbr]
                    #print lower_bound
                    if lower_bound > param.read_len + 20 * k:
#                        for cont_obj in Scaffolds[node[0] ].contigs:
#                            print cont_obj.name
#                        try:
#                            nr_links = G[node][nbr]['nr_links'] 
#                            print 'links: ',nr_links
#                        except KeyError:
#                            try:
#                                nr_links = G[node][nbr]['nr_links']
#                                print 'links: ', nr_links
#                            except KeyError:
#                                print 'a repeat removed'

                        G.remove_edge(node, nbr)
                        pre_filtered += 1



    print 'Nr of edges that did not pass the pre filtering step: ', pre_filtered
    return(G)

def CalculateStats(sorted_contig_lengths, param):
    cur_length = 0
    nr_conts = 0
    LG50 = 0
    NG50 = 0
    for contig_length in sorted_contig_lengths:
        cur_length += contig_length
        nr_conts += 1
        if cur_length >= param.tot_assembly_length / 2.0:
            LG50 = contig_length
            NG50 = nr_conts
            break
    if LG50 == 0:
        print 'NG50 and LG50 could not be calculated in this step data cleared to keep low memory profile.'
        print 'This will be taken care of in later versions.'
    print 'LG50: ', LG50, 'NG50: ', NG50, 'Initial contig assembly length: ', param.tot_assembly_length
    return(NG50, LG50)

def CalculateMeanCoverage(Contigs, first_lib, output_dest, bamfile):
    # tuples like (cont lenght, contig name)
    list_of_cont_tuples = [(Contigs[contig].length, contig) for contig in Contigs]
    #sorted as longest first
    list_of_cont_tuples = sorted(list_of_cont_tuples, key=lambda tuple: tuple[0], reverse=True)
    #coverages of longest contigs
    longest_contigs = list_of_cont_tuples[:1000]
    cov_of_longest_contigs = [Contigs[contig[1]].coverage for contig in longest_contigs]
    #Calculate mean coverage from the 1000 longest contigs
    n = float(len(cov_of_longest_contigs))
    mean_cov = sum(cov_of_longest_contigs) / n
    std_dev = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_cov + mean_cov ** 2), cov_of_longest_contigs))) / (n - 1)) ** 0.5
    extreme_obs_occur = True
    print 'Mean coverage before filtering out extreme observations = ', mean_cov
    print 'Std dev of coverage before filtering out extreme observations= ', std_dev

    ## SMOOTH OUT THE MEAN HERE by removing extreme observations## 
    while extreme_obs_occur:
        extreme_obs_occur, filtered_list = RemoveOutliers(mean_cov, std_dev, cov_of_longest_contigs)
        n = float(len(filtered_list))
        mean_cov = sum(filtered_list) / n
        std_dev = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_cov + mean_cov ** 2), filtered_list))) / (n - 1)) ** 0.5
        cov_of_longest_contigs = filtered_list

    print 'Mean coverage after filtering = ', mean_cov
    print 'Std coverage after filtering = ', std_dev
    print 'Length of longest contig in calc of coverage: ', longest_contigs[0][0]
    print 'Length of shortest contig in calc of coverage: ', longest_contigs[-1][0]


    #fig = plt.figure() 
    #coverage of 1000 longest contigs (unfiltered for extreme observations for the first lib)
    # Why not use filtered list instead!?
    try:
        import matplotlib
        plt.hist(cov_of_longest_contigs, bins=50)
        library = bamfile.split('/')[-1]
        plt.savefig(output_dest + "/BESST_cov_1000_longest_cont" + library + ".png")
        plt.clf()
    except ImportError:
        pass
    return(mean_cov, std_dev)

def RemoveOutliers(mean_cov, std_dev, cov_list):
    k = ChooseQuantileDetectingRepeats(len(cov_list))
    filtered_list = list(filter((lambda x : x < mean_cov + k * std_dev), cov_list))
    if len(cov_list) > len(filtered_list):
        filtered_list = list(filter((lambda x : x < mean_cov + k * std_dev), cov_list))
        return(True, filtered_list)
    else:
        return(False, filtered_list)



def ChooseQuantileDetectingRepeats(nr_of_contigs):
    # Here, we choose the quantile that separates "normal" contigs from repeats. We want to find 
    # the k s.t. the probability of marking one or more normal contigs (from the hole set of contigs) 
    # as repeats (false positives) is less than p=0.05. This is equivalent to: Choosing p in Bin(n,p)
    # (where n = nr of contigs) s.t. P(k=0)>= 0.95 (no successes if a success is a false positive).
    # We get P(k=0)= choose(n,k)*p**n => p = 1 - (0.95*n)**n. With this value of p, if X~N(mean,sigma),
    # we want to derive k from P_x(x < mean + k*sigma) = 1-p. This gives the k that is returned from this function. 
    #from scipy.stats import norm
    import math
    def rational_approximation(t):

        # Abramowitz and Stegun formula 26.2.23.
        # The absolute value of the error should be less than 4.5 e-4.
        c = [2.515517, 0.802853, 0.010328]
        d = [1.432788, 0.189269, 0.001308]
        numerator = (c[2] * t + c[1]) * t + c[0]
        denominator = ((d[2] * t + d[1]) * t + d[0]) * t + 1.0
        return t - numerator / denominator

    def normal_CDF_inverse(p):

        assert p > 0.0 and p < 1

        # See article above for explanation of this section.
        if p < 0.5:
            # F^-1(p) = - G^-1(p)
            return -rational_approximation(math.sqrt(-2.0 * math.log(p)))
        else:
            # F^-1(p) = G^-1(1-p)
            return rational_approximation(math.sqrt(-2.0 * math.log(1.0 - p)))

    p = 1 - (0.95) ** (1 / float(nr_of_contigs))
    #k=norm.ppf(1-p)
    k = normal_CDF_inverse(1 - p)
    print 'Quantile for repeat detector chosen to:', k
    return(k)



def GetParams(bam_file, param, Scaffolds, C_dict, F, Contigs):
    import sys
    informative_pair = set([147, 163]) #161,145,129,177,
    cont_names = bam_file.references
    cont_lengths = bam_file.lengths
    #cont_lengths=[int(nr) for nr in cont_lengths]  #convert long to int object
    cont_lengths_list = list(cont_lengths)
    print "number of contigs:", len(cont_lengths_list)
    indexes = [i for i in range(0, len(cont_lengths_list))]
    from heapq import nlargest
    largest_contigs_indexes = nlargest(1000, indexes, key=lambda i: cont_lengths_list[i]) #get indexes of the 1000 longest contigs
    print largest_contigs_indexes
    if not param.read_len: # user has not specified read len  
        #get read length
        try:
            iter = bam_file.fetch(cont_names[largest_contigs_indexes[0]])
        except ValueError:
            sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
            sys.exit(0)
        nr_reads = 0
        tot_read_len = 0
        for read in iter:
            print read
            if read.rlen != 0:
                tot_read_len += read.rlen
                nr_reads += 1
            else:
                tot_read_len += read.alen
                nr_reads += 1
                #print 'Read has no reported length'
        param.read_len = tot_read_len / float(nr_reads)




    print ''
    print 'Mean of library set to: No mean calc since RNA reads'
    print 'Standard deviation of library set to: No std calc since RNA reads'
    print '-T (library insert size threshold) set to: ', param.ins_size_threshold
    print '-k set to (Scaffolding with contigs larger than): ', param.contig_threshold
    print 'Number of links required to create an edge: ', param.edgesupport
    print 'Read length set to: ', param.read_len
    print 'Relative weight of dominating link set to (default=3): ', param.rel_weight
    print ''
    return()

def PosDirCalculatorPE(cont_dir1, read_dir, cont1pos, readpos, s1len, cont1_len, cont_dir2, mate_dir, cont2pos, matepos, s2len, cont2_len, read_len):
    if cont_dir1 and read_dir:
        obs1 = s1len - cont1pos - readpos
        read_side1 = 'R'
    if cont_dir2 and mate_dir:
        obs2 = s2len - cont2pos - matepos
        read_side2 = 'R'
    if (not cont_dir1) and read_dir:
        obs1 = cont1pos + (cont1_len - readpos)
        read_side1 = 'L'
    if (not cont_dir2) and mate_dir:
        obs2 = cont2pos + (cont2_len - matepos)
        read_side2 = 'L'
    if cont_dir1 and not read_dir:
        obs1 = cont1pos + readpos + read_len
        read_side1 = 'L'
    if cont_dir2 and not mate_dir:
        obs2 = cont2pos + matepos + read_len
        read_side2 = 'L'
    if not cont_dir1 and not read_dir:
        obs1 = s1len - cont1pos - (cont1_len - readpos - read_len)
        read_side1 = 'R'
    if not cont_dir2 and not mate_dir:
        obs2 = s2len - cont2pos - (cont2_len - matepos - read_len)
        read_side2 = 'R'
    #obs=obs1+obs2
    if read_side1 == 'L':
        scaf_side1 = 'L'
    if read_side2 == 'L':
        scaf_side2 = 'L'
    if read_side1 == 'R':
        scaf_side1 = 'R'
    if read_side2 == 'R':
        scaf_side2 = 'R'
    return(obs1, obs2, scaf_side1, scaf_side2)







