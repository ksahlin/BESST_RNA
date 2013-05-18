'''
Created on Sep 23, 2011

@author: ksahlin
'''
import Contig, Scaffold, Parameter
import GenerateOutput as GO
import networkx as nx
from networkx import algorithms
import sys
import GapCalculator as GC
from Norm import normpdf, normcdf


def Algorithm(G, Contigs, Scaffolds, F, Information, C_dict, param):
    #search for linear streches in the graph, remove all cliques >2 and all contigs having more than three neighbors
    nr_edges = 0
    for edge in G.edges_iter():
        if G[edge[0]][edge[1]]['nr_links']:
            nr_edges += 1

    print >> Information, str(nr_edges) + ' link edges created.'
    print 'Perform inference on scaffold graph...'


    #save graph in dot format to file here

    ##If sigma specified. Pre calculate a look up table for every possible gap estimate in the common case
    ##where we have two long contigs that are linked. We do this one time per library to save computation time.
        #if param.std_dev_ins_size:
                #dValuesTable=GC.PreCalcMLvaluesOfdLongContigs(param.mean_ins_size,param.std_dev_ins_size,param.read_len)
        #else :
    dValuesTable = None

    ##### Here is the scaffolding algorithm #######
    ##
    # Pre-scaffolding step here: include small contigs from G' that are occurring in linear regions.
    ##
    G = RemoveIsolatedContigs(G, Information)     #step1
#G = RemoveSpuriousEdges(G,Scaffolds,Contigs,Information,param)    #step2 Remove edges that has very few edges compared to what the coverage and gap estimate says (should I also remove those with too many compared to what GapEst says? -> implement later)
    #G,Contigs,Scaffolds = DetectHaplotypes(G,Contigs,Scaffolds,param,C_dict) #remove haplotypes from scaffolds (to get more linear regions), one haplotype is left in graph and one in outputted to a file haplotypes.fa (with reference to which contig it is haplotypic to) 
    G = RemoveAmbiguousRegions(G, Information, param) #step2
    G = RemoveIsolatedContigs(G, Information) #there are probably new isolated nodes created from step 2
    G, Contigs, Scaffolds = RemoveLoops(G, Scaffolds, Contigs, Information, F)    #step4

    #The contigs that made it to proper scaffolds
    (Contigs, Scaffolds, F, param) = NewContigsScaffolds(G, Contigs, Scaffolds, F, Information, C_dict, dValuesTable, param)   #step5
    ####### End of algorithm #####################

    return(Contigs, Scaffolds, F, param)

def VizualizeGraph(G, param):
    import os
    try:
        import matplotlib
        matplotlib.use('Agg')

        try:
            os.mkdir(param.output_directory + '/graph_regions' + str(int(param.mean_ins_size)))
        except OSError:
            #directory is already created
            pass
        #CB = nx.connected_component_subgraphs(G)
        CB = nx.cycle_basis(G)
        print 'NR of SubG: ' , len(CB)
        counter = 1
        for cycle in CB:
            if len(cycle) >= 6: # at leats 6 nodes if proper cycle (haplotypic region)
                subgraph = nx.Graph()
                subgraph.add_cycle(cycle)
                nx.draw(subgraph)
                matplotlib.pyplot.savefig(param.output_directory + 'graph_regions' + str(int(param.mean_ins_size)) + '/' + str(counter) + '.png')
                matplotlib.pyplot.clf()
                counter += 1
    except ImportError:
        pass
    return()

def RemoveIsolatedContigs(G, Information):
    print 'Remove isolated nodes.'
    counter = 0
    for node in G.nodes():
        if node in G:
            nbr = G.neighbors(node)[0]
            if len(G.neighbors(node)) == 1 and len(G.neighbors(nbr)) == 1:
                counter += 1
                G.remove_nodes_from([node, nbr])
    print >> Information, str(counter) + ' isolated contigs removed from graph.'
    return(G)


def RemoveAmbiguousRegions(G, Information, param):
### remove edges from node if more than two edges 
### Keep the edge with max support e* if all the other edges of that node satisfies links(e*)>l*links(e_i) (l input threshold) all i
### and links(e*)>edge_support (input value)
    edge_support = param.edgesupport
    ratio = param.rel_weight
    print 'Remove edges from node if more than two edges'
    counter1 = 0
    counter2 = 0
    counter3 = 0
    for node in G:
        nbrs = G.neighbors(node)
        #Remove ambiguous edges
        if len(nbrs) > 2:
            nr_link_list = []
            for nbr in nbrs:
                if G[node][nbr]['nr_links']:
                    nr_link_list.append((G[node][nbr]['nr_links'], nbr))
            nr_link_list.sort()

        ###### this if statement is used when scaffolding with fosmidpools, i.e. -p <int> is specified ###
            if param.fosmidpool != None and nr_link_list[-2][0] > param.fosmidpool:
            ### remove all edges on this side of the contig
                nr_nbrs = len(nr_link_list)
                for i in xrange(0, nr_nbrs):
                    G.remove_edge(node, nr_link_list[i][1])
                counter3 += 1
                continue
            ##################################################################

            if nr_link_list[-1][0] >= (ratio * nr_link_list[-2][0]) and nr_link_list[-1][0] >= edge_support:
            ### save the dominating link edge on this side of the contig
                #print 'Linklist: ', nr_link_list
                nr_nbrs = len(nr_link_list)
                for i in xrange(0, nr_nbrs - 1):
                    G.remove_edge(node, nr_link_list[i][1])
                counter1 += 1
            else:
            ### remove all edges on this side of the contig
                nr_nbrs = len(nr_link_list)
                for i in xrange(0, nr_nbrs):
                    G.remove_edge(node, nr_link_list[i][1])
                counter2 += 1

        #Remove low support edges e < k in remaining linear regions             
        else:
            for nbr in nbrs:
                if G[node][nbr]['nr_links'] and G[node][nbr]['nr_links'] < edge_support:
                    G.remove_edge(node, nbr)


    print >> Information, str(counter1 + counter2) + ' ambiguous regions in graph ( a contig with more than 2 neighbors).'
    print >> Information, str(counter1) + ' of these regions solved (one dominating edge)'
    print >> Information, str(counter2) + ' of these regions unsolved (removed all edges)'
    if param.fosmidpool != None:
        print str(counter3), 'of these regions were unsolved due to fosmid pool specific settings (removed all edges)'
    return(G)





def RemoveLoops(G, Scaffolds, Contigs, Information, F):
#### After the proceure above, we hope that the graph is almost perfectly linear but we can still be encountering cycles (because of repeats or haplotypic contigs that has slipped through our conditions). Thus we finally search for loops
    print 'Contigs/scaffolds left:', len(G.nodes()) / 2
    print 'Remove remaining cycles...'
    graphs = nx.connected_component_subgraphs(G)
    #print 'Nr connected components',len(graphs)
    counter = 0
    for graph in graphs:
        list_of_cycles = algorithms.cycles.cycle_basis(graph)
        marked_scaf = set()
        for cycle in list_of_cycles:
            #print >> Information, 'A cycle in the scaffold graph: '+str(cycle)+'\n'
            print >> Information, 'A cycle in the scaffold graph consisting of input contigs/scaffolds: '
            for scaf in cycle:
                if Scaffolds[scaf[0]].contigs[0] not in marked_scaf:
                    print >> Information, Scaffolds[scaf[0]].contigs[0].name
                    marked_scaf.add(Scaffolds[scaf[0]].contigs[0])
            counter += 1
            for node in cycle:
                if node in G:
                    #we split up the whole cycle into separate contigs and send them to F
                    scaffold_ = node[0]
                    G.remove_nodes_from([(scaffold_, 'L'), (scaffold_, 'R')])
                    S_obj = Scaffolds[scaffold_]
                    list_of_contigs = S_obj.contigs   #list of contig objects contained in scaffold object
                    Contigs, F = GO.WriteToF(F, Contigs, list_of_contigs)
                    del Scaffolds[scaffold_]
    print >> Information, str(counter) + ' cycles removed from graph.'
    return(G, Contigs, Scaffolds)

def NewContigsScaffolds(G, Contigs, Scaffolds, F, Information, C_dict, dValuesTable, param):
### Remaining scaffolds are true sensible scaffolds, we must now update both the library of scaffold objects and the library of contig objects
    new_scaffolds_ = nx.connected_component_subgraphs(G)
    print 'Nr of new scaffolds created: ' + str(len(new_scaffolds_))
    print >> Information, 'Nr of new scaffolds created in this step: ' + str(len(new_scaffolds_))
    for new_scaffold_ in new_scaffolds_:
        param.scaffold_indexer += 1
        #scaf_size=len(new_scaffold_)
        scaffold_length = 0
        contig_list = []

        #Store nr_of links between contigs before "destroying" the graph
        for edge in new_scaffold_.edges_iter():
            nr_links = G[edge[0]][edge[1]]['nr_links']
            side1 = edge[0][1]
            side2 = edge[1][1]
            if nr_links:
                contig_objects1 = Scaffolds[edge[0][0]].contigs
                contig_objects2 = Scaffolds[edge[1][0]].contigs
                GiveLinkConnection(Contigs, contig_objects1, contig_objects2, side1, side2, nr_links)

        for node in new_scaffold_:
            if len(G.neighbors(node)) == 1:
                break

        #Create info to new scaffold object such as total length and the contig objects included

        prev_node = ('', '')
        pos = 0
        (G, contig_list, scaffold_length) = UpdateInfo(G, Contigs, Scaffolds, node, prev_node, pos, contig_list, scaffold_length, C_dict, dValuesTable, param)
        S = Scaffold.scaffold(param.scaffold_indexer, contig_list, scaffold_length, {}, {})  #Create the new scaffold object 

        Scaffolds[S.name] = S        #include in scaffold library


    return(Contigs, Scaffolds, F, param)

def UpdateInfo(G, Contigs, Scaffolds, node, prev_node, pos, contig_list, scaffold_length, C_dict, dValuesTable, param):
    scaf = node[0]
    side = node[1]
    prev_scaf = prev_node[0]
    if len(G.neighbors((scaf, side))) == 0:  #reached end of scaffol
        #find the contig with the largest position
        object_with_largest_pos_in_scaffold = max(contig_list, key=lambda object: object.position + object.length)
        scaffold_length = object_with_largest_pos_in_scaffold.position + object_with_largest_pos_in_scaffold.length
        del Scaffolds[scaf] #finally, delete the old scaffold object
        return(G, contig_list, scaffold_length)
    else:
        nbr_node = G.neighbors((scaf, side))
        nbr_scaf = nbr_node[0][0]
        nbr_side = nbr_node[0][1]
        if scaf != prev_scaf:
            if side == 'L':    #Contig/scaffold still has same orientation as in previous iteration, just update position in scaffold                                           
                #want to assign nr of links to contig object, note that in case of a "multiple contigs"-scaffold object, only the outermost contig holds the information of the total nr of links between the two scaffold objects
                contig_objects = Scaffolds[scaf].contigs #list of contig objects

                #Update just update position in scaffold 
                for contig in contig_objects:
                    contig.scaffold = param.scaffold_indexer
                    contig.position += pos
                    #direction unchanged
                    contig_list.append(contig)
                G.remove_node((scaf, side))
                prev_node = node
                node = (nbr_scaf, nbr_side)
                pos += Scaffolds[scaf].s_length  #update position before sending it to next scaffold
                G, contig_list, scaffold_length = UpdateInfo(G, Contigs, Scaffolds, node, prev_node, pos, contig_list, scaffold_length, C_dict, dValuesTable, param)

            else:  #Contig/scaffold need to change orientation as well as modify orientation in this case
                contig_objects = Scaffolds[scaf].contigs #list of contig objects
                for contig in contig_objects:
                    contig.scaffold = param.scaffold_indexer
                    curr_scaf_length = Scaffolds[scaf].s_length
                    curr_pos_within_scaf = contig.position
                    contig_length = contig.length
                    contig.position = pos + (curr_scaf_length - curr_pos_within_scaf) - contig_length #updates the position within scaf
                    contig.direction = bool(True -contig.direction) #changes the direction
                    contig_list.append(contig)

                G.remove_node((scaf, side))
                prev_node = node
                node = (nbr_scaf, nbr_side)
                pos += Scaffolds[scaf].s_length  #update position before sending it to next scaffold
                G, contig_list, scaffold_length = UpdateInfo(G, Contigs, Scaffolds, node, prev_node, pos, contig_list, scaffold_length, C_dict, dValuesTable, param)
        else:
            #calculate gap to next scaffold
            sum_obs = G[(scaf, side)][(nbr_scaf, nbr_side)]['gap_dist']
            nr_links = G[(scaf, side)][(nbr_scaf, nbr_side)]['nr_links']
            #data_observation=(nr_links*param.mean_ins_size -sum_obs)/float(nr_links)
            c1_len = Scaffolds[scaf].s_length
            c2_len = Scaffolds[nbr_scaf].s_length
            avg_gap = 100
            pos += int(avg_gap)
            G.remove_node((scaf, side))
            prev_node = node
            node = (nbr_scaf, nbr_side)
            del Scaffolds[scaf] #finally, delete the old scaffold object
            G, contig_list, scaffold_length = UpdateInfo(G, Contigs, Scaffolds, node, prev_node, pos, contig_list, scaffold_length, C_dict, dValuesTable, param)
    return(G, contig_list, scaffold_length)


def GiveLinkConnection(Contigs, contig_objects1, contig_objects2, side1, side2, nr_links):
    if side1 == 'R' and side2 == 'L':
        max_pos = 0
        for contig in contig_objects1:
            if contig.position >= max_pos:
                linking_contig1 = contig
                max_pos = contig.position
        min_pos = sys.maxint
        for contig in contig_objects2:
            if contig.position <= min_pos:
                linking_contig2 = contig
                min_pos = contig.position

        linking_contig1.links[linking_contig2.name] = nr_links
        linking_contig2.links[linking_contig1.name] = nr_links

        #print linking_contig1.name, linking_contig2.name, nr_links


    elif side1 == 'L' and side2 == 'R':
        max_pos = 0
        for contig in contig_objects2:
            if contig.position >= max_pos:
                linking_contig2 = contig
                max_pos = contig.position

        min_pos = sys.maxint
        for contig in contig_objects1:
            if contig.position <= min_pos:
                linking_contig1 = contig
                min_pos = contig.position
        linking_contig1.links[linking_contig2.name] = nr_links
        linking_contig2.links[linking_contig1.name] = nr_links

        #print linking_contig1.name, linking_contig2.name,nr_links

    elif side1 == 'R' and side2 == 'R':
        max_pos = 0
        for contig in contig_objects1:
            if contig.position >= max_pos:
                linking_contig1 = contig
                max_pos = contig.position

        max_pos = 0
        for contig in contig_objects2:
            if contig.position >= max_pos:
                linking_contig2 = contig
                max_pos = contig.position
        linking_contig1.links[linking_contig2.name] = nr_links
        linking_contig2.links[linking_contig1.name] = nr_links

        #print linking_contig1.name, linking_contig2.name , nr_links

    elif side1 == 'L' and side2 == 'L':
        min_pos = sys.maxint
        for contig in contig_objects1:
            if contig.position <= min_pos:
                linking_contig1 = contig
                min_pos = contig.position
        min_pos = sys.maxint
        for contig in contig_objects2:
            if contig.position <= min_pos:
                linking_contig2 = contig
                min_pos = contig.position

        linking_contig1.links[linking_contig2.name] = nr_links
        linking_contig2.links[linking_contig1.name] = nr_links


        #print linking_contig1.name, linking_contig2.name , nr_links
