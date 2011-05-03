# Input: Connected graph "G", parameter "all"
# Output: A maximal induced tree or the list of all maximal induced trees

from sage.all import *

def max_induced_tree(G,all=False):
    # Return a maximal induced tree, useful for computing tree size of a graph
    if all == False:
        if G.is_tree() == True:
            return G
        else:
            F = Set(G.vertices())
            T = False
            i = 1
            while i < G.num_verts() and T == False:
                card = (F.subsets(i)).cardinality()
                V = list(F.subsets(i))
                j = 0
                while j < card and T==False:
                    G1 = G.copy()
                    G1.delete_vertices(list(V[j]))
                    if G1.is_tree() == True:
                        T = True
                    j = j+1
                i = i+1
            return G1
                   
    else:  # Return list of all maximal induced trees 
        trees = list([])
        if G.is_tree() == True:
            trees.append(G)
            return trees
        else:
            F = Set(G.vertices())
            T = False
            i = 1
            while i < G.size() and T == False:
                card = (F.subsets(i)).cardinality()
                V = list(F.subsets(i))
                j = 0
                while j < card:
                    G1 = G.copy()
                    G1.delete_vertices(list(V[j]))
                    if G1.is_tree() == True:
                        trees.append(G1)
                        T = True
                    j = j+1
                i = i+1
            return trees
			
			
			
def conn_edges_path_verts(G,v1,v2):
    paths = G.all_paths(v1,v2)
    Edges = []
    for i in range(0,len(paths)):
        apath = paths[i]
        for j in range(0, len(apath)-1):
            a = apath[j]
            b = apath[j+1]
            edge = (min(a,b),max(a,b))
            Edges.append(edge)
    return Set(Edges)

# Determine if the following condition holds for a graph G
# There exists a maximal induced tree T of G such that for all v and w not in T, epsilon(v) and epsilon(w) have 
# a nonempty intersection if and only if v and w are adjacent in G

# current code

def msr_is_ts(G, noTree=True):
    
    trees = max_induced_tree(G,all)
    ts = trees[0].num_verts()
    emptySet = (Set([1]).subsets()).first()
    treeFound = False
    goAhead = True
    
    if G.num_verts() - ts < 2:
        goAhead = False
    
    n = 0
    while n < len(trees) and treeFound == False and goAhead == True:
    
        T = trees[n]
        verts = list(Set(G.vertices()).difference(Set(T.vertices()))) # Vertices of G not in T
        ext_pairs = list((Set(verts)).subsets(2)) # Iterate over all pairs of vertices in verts
        keepGoing = 1 # If given any tree, epsilon(v) cap epsilon(w) = {} iff v,w adjacent fails, move to next tree
        
        k = 0
        while k < len(ext_pairs) and keepGoing == 1:
            
            twoVerts = list([list(ext_pairs[k])[0],list(ext_pairs[k])[1]]) #put v,w in list form
            twoVertsEdgeSets = list([])
           
            for i in range(0,2): # find adjacent vertices of v and w in T
            
                adjacent = list([])
                for j in T.vertices():
                    if (min(twoVerts[i],j),max(twoVerts[i],j),{}) in G.edges():
                        adjacent.append(j)
                adjacent = Set(adjacent)
                int_pairs = list(adjacent.subsets(2)) # Iterate over all vertices in T adjacent to v (or w)
                    
                for l in range(0,len(int_pairs)): #build epsilon(v) and epsilon(w)
                
                    v1 = list(int_pairs[l])[0]
                    v2 = list(int_pairs[l])[1]
                    if l == 0:
                        edgeSets = conn_edges_path_verts(T, v1,v2)
                    else:
                        edgeSets = edgeSets.union(conn_edges_path_verts(T,v1,v2))
                        
                if len(int_pairs) == 0:
                    twoVertsEdgeSets.append(emptySet)
                else:
                    twoVertsEdgeSets.append(edgeSets)
                    
            if (min(twoVerts[0],twoVerts[1]),max(twoVerts[0],twoVerts[1]),{}) in G.edges(): # v,w adjacent
                
                if twoVertsEdgeSets[0].intersection(twoVertsEdgeSets[1]) == emptySet: # iff fails
                    keepGoing = 0
                         
            else: # v,w not adjacent
                
                if twoVertsEdgeSets[0].intersection(twoVertsEdgeSets[1]) != emptySet: # iff fails
                    keepGoing = 0  
                        
            k = k + 1
            
        n = n + 1    
        
        if keepGoing == 1:
            treeFound = True
            success = 'YES'
           
        if n == len(trees) and treeFound == False:
            success = 'NO'
    
    if goAhead == False:
        print 'YES'
    else:
        print success
