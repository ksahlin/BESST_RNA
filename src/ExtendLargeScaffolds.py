'''
Created on Jun 4, 2012

@author: ksahlin
'''

from time import time
import networkx as nx

G = nx.Graph() 
#G.add_nodes_from([(1, 'L'), (1, 'R'), (2, 'L'), (2, 'R'), (3, 'L'), (3, 'R'), (4, 'L'), (4, 'R'), (5, 'L'), (5, 'R')]) 
for i in range(1, 7):
    G.add_edge((i, 'L'), (i, 'R'), {'weight':0})
G.add_edges_from([((1, 'R'), (2, 'R'), {'weight':1}), ((3, 'L'), (4, 'L'), {'weight':1}), ((2, 'L'), (3, 'R'), {'weight':1}), ((1, 'R'), (5, 'L'), {'weight':2}),
                   ((5, 'R'), (4, 'L'), {'weight':3}), ((2, 'L'), (5, 'L'), {'weight':2}), ((1, 'R'), (4, 'L'), {'weight':8}), ((2, 'L'), (6, 'L'), {'weight':3}),
                   ((1, 'L'), (4, 'R'), {'weight':1}), ((1, 'L'), (4, 'L'), {'weight':1}), ((3, 'L'), (4, 'R'), {'weight':1})]) 
G_prime = nx.Graph()
G_prime.add_nodes_from([(1, 'L'), (1, 'R'), (4, 'L'), (4, 'R'), (6, 'L'), (6, 'R')])
contigs = [1, 2, 3, 4, 5, 6]


def find_all_paths_between_fixed_nodes(graph, start, end): 
    path = [] 
    paths = [] 
    if start[1] == 'L' and end[1] == 'L':
        forbidden = [(start[0], 'R'), (end[0], 'R')]
    elif start[1] == 'L' and end[1] == 'R':
        forbidden = [(start[0], 'R'), (end[0], 'L')]
    elif start[1] == 'R' and end[1] == 'L':
        forbidden = [(start[0], 'L'), (end[0], 'R')]
    else:        
        forbidden = [(start[0], 'L'), (end[0], 'L')]
    sum_path = 0
    queue = [(start, end, path, sum_path)]
    prev_node = start
    while queue: 
        start, end, path, sum_path = queue.pop()            
        #print 'PATH', path 
        path = path + [start] 
        if sum_path > 19: #replace with shortest path
            continue
        if start == end:
            print 'PATH', path, sum_path
            paths.append(path) 
            continue          
        if  prev_node[0] != start[0]:
            if start[1] == 'L' and (start[0], 'R') not in forbidden: 
                queue.append(((start[0], 'R'), end, path, sum_path + graph[start][(start[0], 'R')]['weight']))
            elif start[1] == 'R' and (start[0], 'L') not in forbidden: 
                queue.append(((start[0], 'L'), end, path, sum_path + graph[start][(start[0], 'L')]['weight']))                
        else:
            for node in set(graph[start]).difference(path): 
                if node not in forbidden: 
                    queue.append((node, end, path, sum_path + graph[start][node]['weight'])) 
        prev_node = start        
    return paths

def Insert_path(a, x, path, lo=0, hi=None): 
    if len(a) == 0:
        a.append([x, path])
        return()        
    if x > a[-1][0]:
        a.append([x, path])
        return()
    if x < a[0][0]:
        a.insert(0, [x, path])
        return()
    if hi is None: 
        hi = len(a) 
    while lo < hi: 
        mid = (lo + hi) // 2 
        midval = a[mid][0] 
        if midval < x: 
            lo = mid + 1 
        elif midval > x:  
            hi = mid 
        else: 
            #exact same score as some other path, look what path has the most nr of contigs and prefer that one
            if len(a[mid][1]) > len(path):
                a.insert(mid, [x, path])
                print mid
            else:
                a.insert(mid + 1, [x, path])
            return () 
    if x > a[lo]:
        a.insert(lo + 1, [x, path])
        print mid
    else:
        a.insert(lo, [x, path])
        print mid
    return () 

def ScorePaths(G, nodes_present_in_path, paths, all_paths_sorted_wrt_score):
    if len(paths) == 0:
        return 0
    def CalculateSpanScore(path, G): 
        contig_end_prev = (-1, 'O')
        contig_start_prev = (-1, 'O') 
        spanned_contigs_in_path = set([])
        for j, contig_end in reversed(list(enumerate(path))):
            if contig_end_prev[0] == contig_end[0]:
                continue
            for i, contig_start in reversed(list(enumerate(path))):
                if contig_start_prev[0] == contig_start[0]:
                    continue   
                if contig_end in G[contig_start]:
                    for contig in path[i + 1:j]:
                        spanned_contigs_in_path.add(contig[0])
                           
            contig_end_prev = contig_end
            contig_start_prev = contig_start
            
        #print spanned_contigs_in_path
        score = len(spanned_contigs_in_path)        
        return(score)
    
    def CalculateNbrScore(path, nodes_present_in_path, path_start, path_end, G):
        tot_nbr_ratio = 0
        #print path
        tot_nbrs = 0
        tot_nbrs_not_in_path = 0
        edges_already_considered = set()
        for contig_end in path: #[1:-1]:
            for nbr in G.neighbors(contig_end):
                if nbr[0] != contig_end[0] and (nbr, contig_end) not in edges_already_considered:
                    tot_nbrs += 1
                    edges_already_considered.add((contig_end, nbr))
                    edges_already_considered.add((nbr, contig_end))
                    if nbr not in nodes_present_in_path[(path_start, path_end)]:
                        tot_nbrs_not_in_path += 1                    
        if tot_nbrs != 0:  
            #print tot_nbrs, tot_nbrs_not_in_path
            tot_nbr_ratio += (tot_nbrs - tot_nbrs_not_in_path) / float(tot_nbrs)    
        else:
            #direct path between two scaffolds with larger contigs
            print 'These contigs/scaffolds should have been joined earlier.. why did we get here?'
            
        return tot_nbr_ratio
    
    for path in paths:
        path_start = path[0]
        path_end = path[-1]
        #calculate spanning score s_ci
        span_score = CalculateSpanScore(path, G)
        #calculate neighbour explaining path score r_ci
        nbr_score = CalculateNbrScore(path, nodes_present_in_path, path_start, path_end, G)
        #print nbr_score, span_score
        if len(path) > 2: #startnode and end node are not directly connected
            tot_score = (span_score) / ((len(path) - 2) / 2.0) + nbr_score
        else:
            tot_score = 0 #if they are directly connected, they should by definition already be in scaffold together
            
        #insert path in sorted list using binary search
        Insert_path(all_paths_sorted_wrt_score, tot_score, path)

    return ()

def find_all_paths_for_start_node(graph, start, end, nodes_present_in_path, already_visited): 
    path = [] 
    paths = [] 
    if start[1] == 'L':
        forbidden = (start[0], 'R')
    else:
        forbidden = (start[0], 'L')
    start_node = start
#    for end_node in end.difference(already_visited.union([forbidden])):
#        nodes_present_in_path[(start_node, end_node)] = set()
    sum_path = 0
    queue = [(start, end, path, sum_path)]
    prev_node = start
    while queue: 
        start, end, path, sum_path = queue.pop()            
        #print 'PATH', path 
        path = path + [start] 
        #if sum_path > 19: #All possible paths can be exponential!! need something to stop algorithm in time
        #    continue
        #if score < score_best_path: # need something to stop a bad path
        #    continue
        if start in already_visited or start in forbidden:
            continue
        if start in end:
            if (start_node, start) in nodes_present_in_path:
                nodes_present_in_path[(start_node, start)] = nodes_present_in_path[(start_node, start)].union(path)
            else:
                nodes_present_in_path[(start_node, start)] = set(path)
                
            #print 'PATH', path, sum_path
            #score = ScorePath(G, path)
            paths.append(path) 
            #nodes_present_in_path[(start_node, start)] = nodes_present_in_path[(start_node, start)].union(path)
            continue          
        if  prev_node[0] != start[0]:
            if start[1] == 'L' and (start[0], 'R') != forbidden: 
                queue.append(((start[0], 'R'), end, path, sum_path + graph[start][(start[0], 'R')]['weight']))
            elif start[1] == 'R' and (start[0], 'L') != forbidden: 
                queue.append(((start[0], 'L'), end, path, sum_path + graph[start][(start[0], 'L')]['weight']))                
        else:
            for node in set(graph[start]).difference(path): 
                if node != forbidden: 
                    queue.append((node, end, path, sum_path + graph[start][node]['weight'])) 
        prev_node = start        
    return paths

def ExtendScaffolds(all_paths_sorted_wrt_score):
    for score_and_path in reversed(all_paths_sorted_wrt_score):
        print 'Score: ',score_and_path[0],'Path: ' , score_and_path[1]
        #MakeScaffolds()?
    return()


start = time()
paths = find_all_paths_between_fixed_nodes(G, (1, 'R'), (4, 'L')) 
elapsed = time() - start
print len(paths)
for path in paths:
    print path
print 'time all paths: ', elapsed


def BetweenScaffolds(G,G_prime):
    start = time()
    end = set()
    for node in G_prime:
        end.add(node)
    
    # here we should have a for loop looping over all start nodes. Start nodes already examined should be removed in a nice way to skip over counting
    iter_nodes = end.copy()
    already_visited = set()
    
    all_paths_sorted_wrt_score = []
    while len(iter_nodes) > 1: 
        start_node = iter_nodes.pop()
        #print 'START NODE: ', start_node 
        nodes_present_in_path = {}
        paths = find_all_paths_for_start_node(G, start_node, end.difference(set([start_node])), nodes_present_in_path, already_visited)
        already_visited.add(start_node) 
        elapsed = time() - start
        #print 'START NODE: ', start_node, 'Tot nr of paths for this start node: ', len(paths)
        ScorePaths(G, nodes_present_in_path, paths, all_paths_sorted_wrt_score)
        #print  all_paths_sorted_wrt_score
    ExtendScaffolds(all_paths_sorted_wrt_score)
    return()
    

BetweenScaffolds(G,G_prime) 
print 'time all paths: ', elapsed


