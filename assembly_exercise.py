from collections import Counter
import copy

from PRG.exercise_regex.regex_exercise import results


def build_adjacency_list(edges, directed=False):
    """
    Creates the adjacency list of directed or undirected graph
    :param edges: list of 2-tuples representing edges
    :param directed: boolean, True for directed graphs, False for undirected
    :return: adj_list, dict[str, list[str]]
    """
    vertex = []
    for edge in edges:
        vertex.append(edge[0])
        vertex.append(edge[1])
    vertex_unique = list(dict.fromkeys(vertex))
    adj_list = {i: [] for i in vertex_unique}
    for edge in edges:
        u = edge[0]
        v = edge[1]
        adj_list[u].append(v)
        if not directed:
            adj_list[v].append(u)
    return adj_list


def degrees(graph):
    """
    Computes the degree of every vertex of an undirected graph
    :param graph: adjacency list of an undirected graph
    :return: dictionary of vertices and their degrees
    """
    degree_dict = {}
    for vertex in graph:
        degree_dict[vertex] = len(graph[vertex])
    return degree_dict


def in_out_degree(graph):
    """
    Computes the in and out degrees of every vertex of a directed graph
    :param graph: adjacency list of a directed graph
    :return: (in_degree, out_degree), two dict[str, int], keyed by every vertex.
    """
    in_degree = {}
    out_degree = {}
    for vertex in graph:
        out_degree[vertex] = len(graph[vertex])
        in_degree[vertex] = 0
        for vtx in graph:
            in_degree[vertex] = in_degree[vertex] + graph[vtx].count(vertex)
    return in_degree, out_degree


def dfs(graph, start):
    """
    Finds all vertices reachable from the start
    :param graph: adjacency list od a graph (undirected or directed)
    :param start: the starting vertex
    :return: list[str], the vertices reachable from start, in the order you first visit them.
    """
    actual_node = start
    visited = [start] # a set of all visited vertices
    walk = [start] # saving every step I take, so that I can determine, where to get back, when stuck
    while True:
        unvisited_neighbours = [neighbour for neighbour in graph[actual_node] if neighbour not in visited]
        if not unvisited_neighbours: # the case, when I'm stuck
            if actual_node == start: # end the search
                break
            else: # return to the node I initially came here from (using walk)
                idx = walk.index(actual_node) # return the index of the first occurrence of actual_node in the walk
                actual_node = walk[idx-1]
                walk.append(actual_node)
        else:
            unvisited_neighbours.sort()
            actual_node = unvisited_neighbours[0] # move to the alphabetically first unvisited neighbour
            walk.append(actual_node)
            visited.append(actual_node)
    return visited


def connected_components(graph):
    connected_comp =[]
    remaining_vertices = set(graph.keys()) # the vertices, that are not yet in any connected component
    while remaining_vertices:
        start = remaining_vertices.pop()
        component = set(dfs(graph, start))
        connected_comp.append(component) # add the connected component to the list of all connected components
        remaining_vertices = remaining_vertices - component
    return connected_comp


def is_weakly_connected(graph):
    shadow_graph = copy.deepcopy(graph)
    for vertex in graph:
        for destination in graph[vertex]:
            shadow_graph[destination].append(vertex)
    components = connected_components(shadow_graph)
    if len(components) == 1:
        return True
    else:
        return False


def eulerian_type(graph, directed=False):
    result = 'none'
    if directed: # directed
        if is_weakly_connected(graph):
            in_deg, out_deg = in_out_degree(graph)
            deg_plus = 0
            deg_minus = 0
            for node in in_deg:
                if in_deg[node] == out_deg[node]:
                    continue
                elif in_deg[node] == out_deg[node] + 1:
                    deg_plus += 1
                    continue
                elif in_deg[node] == out_deg[node] - 1:
                    deg_minus += 1
                    continue
                else:
                    return result
            if (deg_minus == 0) and (deg_plus == 0):
                result = 'circuit'
            elif (deg_plus == 1) and (deg_minus == 1):
                result = 'path'
            # print(in_deg, out_deg)
    else: # undirected
        if len(connected_components(graph)) == 1:
            degrees_dict = degrees(graph)
            odd_deg = 0
            for node in degrees_dict:
                if (degrees_dict[node] % 2) == 1:
                    odd_deg += 1
            if odd_deg == 0:
                result = 'circuit'
            elif odd_deg == 2:
                result = 'path'
            # print(degrees_dict, odd_deg)
    return result


def kmer_composition(text, k):
    composition = []
    if (k > len(text)) or (k < 1):
        print('Not valid k')
    else:
        for idx in range(len(text)-k+1):
            composition.append(text[idx:idx+k])
    return composition


def debruijn_graph(kmers):
    # divides k-mers into (k-1)-mers only
    result_graph = {}
    k = len(kmers[0])
    for kmer in kmers:
        composition = kmer_composition(kmer, k-1)
        for comp in composition:
            if comp not in result_graph:
                result_graph[comp] = []
        result_graph[composition[0]].append(composition[1])
    return result_graph


def spell_path(nodes):
    path = ''
    if not nodes:
        print('Enter a non-empty list of nodes')
    else:
        for node in nodes:
            if path == '':
                path = node
            else:
                path = path + node[-1]
    return path


def hierholzer(graph, start, verbose=False):
    unused_edges = copy.deepcopy(graph)
    for edge in unused_edges:
        unused_edges[edge].sort(reverse=True) # sort decreasingly
    final_walk = [] # the final walk through all edges
    node = start # the node I am in at the moment
    walk = [node]  # the active walk
    all_walks = [] # list of all walks
    all_splices = []
    while True:
        if unused_edges[node]: # there are unused edges from this node
            node = unused_edges[node].pop() # remove and return the lask index, tj. alphabeticaly the first
            walk.append(node)
            # mark the edge used
        else: # stuck
            all_walks.append(walk)
            if final_walk: # replace the first occurance of node by the walk
                i = final_walk.index(node)
                final_walk[i:i+1] = walk
            else:
                final_walk = walk
            all_splices.append(final_walk.copy())
            if not any(unused_edges.values()): # the algorithm is finished
                break
            for vertex in final_walk:
                if unused_edges[vertex]:
                    node = vertex
                    break
            walk = [node]

    if verbose:
        return final_walk, all_walks, all_splices
    else:
        return final_walk


def eulerian_path(graph, verbose=False):
    # find the starting node
    start = ''
    graph_type = eulerian_type(graph, directed=True)
    if graph_type == 'circuit':
        start = min(graph)
    elif graph_type == 'path':
        in_deg, out_deg = in_out_degree(graph)
        for node in graph:
            if out_deg[node] > in_deg[node]:
                start = node
    else:
        print('The Eulerian path nor circuit do not exist.')
        return 1

    # get the eulerian path
    if verbose:
        path, all_paths, all_splices = hierholzer(graph, start, verbose)
        result = (spell_path(path), all_paths, all_splices)
    else:
        path = hierholzer(graph, start, verbose)
        result = spell_path(path)

    return result

if __name__ == "__main__":
    de_graph = debruijn_graph(['ACT', 'AGA', 'CAG', 'CCA', 'CTG', 'GAA', 'GAC', 'GGT', 'GTC', 'TCC', 'TGA'])
    print(de_graph)
    print(spell_path(["GAT", "ATT", "TTA"]))
    print(spell_path(['GG', 'GT', 'TC', 'CC', 'CA', 'AG', 'GA', 'AC', 'CT', 'TG', 'GA', 'AA']))

    # Hierholzer's algotihm
    de_graph1 = debruijn_graph(['ACT', 'AGC', 'ATT', 'CCG', 'CGA', 'CGT', 'CTG', 'GAC', 'GAG', 'GAT', 'GTC', 'TCC', 'TCG', 'TGA', 'TTC'])
    print(eulerian_type(de_graph1, directed=True))
    print(hierholzer(de_graph1, 'GA', verbose=True))
    print(eulerian_path(de_graph1), len(eulerian_path(de_graph1)))