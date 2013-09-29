#
#  msr_program.py
#  
#
#  Created by Steven Osborne on 5/25/11.
#  Revised on 7/11/2011.
#  Copyright (c) 2011 Iowa State University. All rights reserved.
#
#
# The source code for Zplus, edge_clique_cover, and find_cut_vertex was authored by 
# (in alphabetical order): Steve Butler, Laura DeLoss, 
# Jason Grout, Tracy Hall, Josh Lagrange, Tracy McKay, Jason Smith, and Geoff Tims. &nbsp;
# It may be found on Github, https://github.com/jasongrout/minimum_rank# 
# This code is distributed under the GNU Public License (GPL), version 2 or (at your option) 
# any later version. See "http://www.gnu.org/licenses/gpl-2.0.html for details

def Zplus_gen(G):
	"""
    Return the positive semidefinite zero forcing number and a positive 
	semidefinite zero forcing set of a given simple graph G.

    INPUT:
        G	--	the graph

    OUTPUT:
        Set containing the positive semidefinite zero forcing number and a positive 
		semidefinite zero forcing set.

    EXAMPLES:
		sage: Zplus_gen(atlas_graphs[0])
		(0, set([]))
        sage: Zplus_gen(graphs.PathGraph(5))
		(1, set([0]))
        sage: Zplus_gen(graphs.CompleteGraph(5))
        (4, set([0, 1, 2, 3]))
        sage: Zplus_gen(graphs.HexahedralGraph())
        (4, set([4, 5, 6, 7]))
    """
	"""
    if G.is_connected()==True:
        if G.num_verts() == 0:
            Zp = 0
            ZFS = []
            return Zp,set(ZFS)
		#All trees have Z+ of 1
        if G.is_tree()==True:
            Zp = 1
            ZFS = [G.vertices()[0]]
            return Zp,set(ZFS)
		#The parent code for Z+ can handle non-tree connected graphs
        else:
            return Zplus(G)
			
    #If G is disconnected, compute Z+ for all components
    else: 
        B = G.connected_components()
        Zp=0
        ZFS=set([])
        #Form each component
        for i in range(0,G.connected_components_number()):
            Geye=G.copy()
            #Delete vertices from other components
            for k in range(0,G.connected_components_number()):
                if k != i:
                    Geye.delete_vertices(B[k])
            #Sum over each component
            if Geye.is_tree() == True:
                S = Geye.vertices()[0]
                Zp = Zp + 1
                ZFS = ZFS.union(set([S]))
            else:
                (Zm,S) = Zplus(Geye)
                Zp = Zp + Zm
                ZFS = ZFS.union(S)
        return Zp,ZFS
	"""
	if G.num_verts() == 0:
		return [0]
	if G.num_verts() == 1:
		return [1]
	else:
		return [Zq_compute(G,0)]

def find_cut_vertex(graph):
    """
    Return a "good" cut-vertex for a graph if it exists; otherwise,
    returns False.

    INPUT:
        graph -- the graph on which to find a cut-vertex

    OUTPUT:
        a cut-vertex (if one exists) that either results in components
        of order less than 7 or a minimum of the maximum component
        order; otherwise False

    EXAMPLES:
        sage: find_cut_vertex(graphs.PathGraph(3))
        1
        sage: find_cut_vertex(graphs.PathGraph(20))
        9
        sage: [find_cut_vertex(graphs.PathGraph(i)) for i in [1..20]]
        [False, False, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 6, 7, 7, 8, 8, 9, 9]
        sage: find_cut_vertex(graphs.CompleteGraph(3))
        False
    """

    vertices=graph.vertices()
    graph_cc_num=graph.connected_components_number()
    graph_order=graph.order()
    
    #this will hold the "best" cut-vertex and the order of the largest
    #connected component after deletion
    best_v=(False,graph_order)

    #checks each vertex and determines the best one
    for v in vertices:
        g=graph.copy()
        g.delete_vertex(v)
        g_cc = g.connected_components()
        if len(g_cc)>graph_cc_num:
            # We have a cut-vertex
            max_order = max(len(c) for c in g_cc)

            if max_order<7:
                return v
            if max_order<best_v[1]:
                best_v=(v,max_order)
            
    return best_v[0]

def max_induced_tree(G,all=False):
	"""
    Return a maximal induced tree, useful for computing tree size of a graph.

    INPUT:
        G	--	the graph
        all --	request to return all maximal induced trees 
				(defualt set to False)

    OUTPUT:
        If all is False, return a maximal induced tree.
		        
       Else, return a list of all the maximal induced trees.

    EXAMPLES:
		sage: max_induced_tree(graphs.PathGraph(3))
		Path Graph: Graph on 3 vertices
		max_induced_tree(graphs.CompleteGraph(5))
		Complete graph: Graph on 2 vertices
		sage: max_induced_tree(atlas_graphs[931])
		Graph on 4 vertices
	"""
	# Return only one maximal induced tree
    if all == False:
		# Take care of trivial case
        if G.is_tree() == True:
            return G
		# Starting with subsets of the vertex set of the graph 
		# of size |G|-1, iterate over successively smaller 
		# vertex sets until a given cardinality admits an 
		# induced tree
        else:
            F = Set(G.vertices())
            T = False
            i = 1
            while i < G.num_verts() and T == False:
                card = (F.subsets(i)).cardinality()
                V = list(F.subsets(i))
                j = 0
				# Test to see if induced sugraph is a tree
                while j < card and T==False:
                    G1 = G.copy()
                    G1.delete_vertices(list(V[j]))
                    if G1.is_tree() == True:
                        T = True
                    j = j+1
                i = i+1
            return G1
			
	# Return list of all maximal induced trees 
    else:
        trees = list([])
        if G.is_tree() == True:
            trees.append(G)
            return trees
        else:
            F = Set(G.vertices())
            T = False
            i = 1
			# Form the list of maximal induced trees
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
	"""
	Return the edge set of all paths between two 
	vertices v1 and v2 of a given graph G.
	
    INPUT:
        G		--	the graph
		v1,v2	--	the vertices in question

    OUTPUT:
       If both v1 and v2 are distinct vertices in the graph G, 
	   return set of edges of the graph G which compose all 
	   possible paths from v1 to v2 in G.
	   Else, return 'At least one of the given vertices
	   is not in the given graph' or 'Provide two distinct 
	   vertices' depending on the circumstance.

    EXAMPLES:
	sage: conn_edges_path_verts(graphs.CompleteGraph(5),1,2)
	{(0, 1), (1, 2), (1, 3), (1, 4), (0, 2), (2, 3), 
	(0, 4), (0, 3), (3, 4), (2, 4)}
	sage: conn_edges_path_verts(atlas_graphs[15],1,4)
	'At least one of the given vertices is not in the given graph'
	"""
	if v1 in G.vertices() and v2 in G.vertices() and v1 != v2:
		# Retrieve a list of the vertex sets that form all the 
		# paths between v1 and v2 in G 
		paths = G.all_paths(v1,v2)
		Edges = []
		# Form edge set admitted by 'paths'
		for i in range(0,len(paths)):
			apath = paths[i]
			for j in range(0, len(apath)-1):
				a = apath[j]
				b = apath[j+1]
				edge = (min(a,b),max(a,b))
				Edges.append(edge)
		return Set(Edges)
	# Error messages
	else: 
		if v1 == v2:
			return 'Provide two distinct vertices'
		else:
			return 'At least one of the given vertices is not in the given graph'

def msr_is_ts(H, noTree=True):
	"""
	Determine if the following condition holds for a graph H:
	There exists a maximal induced tree T of H such that for all 
	v and w not in T, epsilon(v) and epsilon(w) have a nonempty 
	intersection if and only if v and w are adjacent in H.
	
    INPUT:
        H		--	the graph
		noTree	--	request to NOT return the tree that meets the 
					above conditions (default is True)

    OUTPUT:
		If noTree is True and condition holds, return string 'YES', 
			and tree size of the graph (i.e. size of maximal induced tree).
		If noTree is False and condition holds, return string 'YES', 
			tree size of the graph, and a maximal induced tree
			for which the condition holds.
		If the condition fails, return string 'NO' and None.
		
    EXAMPLES:
		sage: msr_is_ts(graphs.CompleteGraph(5))
		('YES', 2)
		sage: msr_is_ts(atlas_graphs[574])
		('NO', None)
		sage: msr_is_ts(Graph({1:[2,3,6],2:[3,4,7],
			3:[4,6,7],4:[5],5:[6,7],6:[7]}),False)
		('YES', 4, Graph on 4 vertices)
	"""
    
    G = Graph(H.adjacency_matrix())
    trees = max_induced_tree(G,all)
    emptySet = (Set([1]).subsets()).first()
    ts = trees[0].num_verts()
    
    treeFound = False
    goAhead = True
    
    if G.num_verts() - ts < 2:
        goAhead = False
    
    n = 0
    while n < len(trees) and treeFound == False and goAhead == True:
    
        T = trees[n]
		# Vertices of G not in T
        verts = list(Set(G.vertices()).difference(Set(T.vertices())))
		 # Iterate over all pairs of vertices in verts
        ext_pairs = list((Set(verts)).subsets(2))
		# If given any tree, epsilon(v) cap epsilon(w) = {} iff v,w adjacent fails, move to next tree
        keepGoing = 1         
        k = 0
        while k < len(ext_pairs) and keepGoing == 1:
            
            twoVerts = list([list(ext_pairs[k])[0],list(ext_pairs[k])[1]]) #put v,w in list form
            twoVertsEdgeSets = list([])
			
            # find adjacent vertices of v and w in T
            for i in range(0,2):
            
                adjacent = list([])
                for j in T.vertices():
                    if (min(twoVerts[i],j),max(twoVerts[i],j),None) in G.edges():
                        adjacent.append(j)
                adjacent = Set(adjacent)
				# Iterate over all vertices in T adjacent to v (or w)
                int_pairs = list(adjacent.subsets(2)) 
                
				# build epsilon(v) and epsilon(w)
                for l in range(0,len(int_pairs)): 
                
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
            
			# v,w adjacent
            if (min(twoVerts[0],twoVerts[1]),max(twoVerts[0],twoVerts[1]), None) in G.edges(): 
                # iff fails
                if twoVertsEdgeSets[0].intersection(twoVertsEdgeSets[1]) == emptySet: 
                    keepGoing = 0
            # v,w not adjacent             
            else: 
                # iff fails
                if twoVertsEdgeSets[0].intersection(twoVertsEdgeSets[1]) != emptySet: 
                    keepGoing = 0  
                        
            k = k + 1
            
        n = n + 1    
        
        if keepGoing == 1:
            treeFound = True
			if noTree == True:
				success = ('YES',ts)
			else: 
				success = ('YES',ts,T)
			
        if n == len(trees) and treeFound == False:
            success = ('NO',None)
    
    if goAhead == False:
		if noTree == True:
			return ('YES',ts)
		else:
			return ('YES',ts,trees[0])
			
    else:
        return success

def cut_vertex_machine(H):
	"""
	Runs the cut-vertex algorithm for msr_bounds.

    INPUT:
        H			--	the graph

    OUTPUT:
       If H is complete or kappa(H) is not 1, 
	   return 'No cut vertex'.
	   Else, return lowerbound for M+(H) using the 
	   cut-vertex algorithm in conjunction with msr_bounds
	   and return explanation.

    EXAMPLES:
	sage: cut_vertex_machine(graphs.PathGraph(5))
	(1, 'Z+(G) <= 3, Z+(G) <= 3, ')
	sage: cut_vertex_machine(graphs.CompleteGraph(5))
	'No cut vertex'
	"""

    G = Graph(H.adjacency_matrix())
    K = graphs.CompleteGraph(G.num_verts())
    if G.is_isomorphic(K) == False:
        kappa = int(G.vertex_connectivity())
        if kappa == 1:
            cut_v = find_cut_vertex(G)
            H = [G.copy(),G.copy(),G.copy()]
            H[0].delete_vertex(cut_v)
            H[1].delete_vertices(H[0].connected_components()[0])
            H[2].delete_vertices(H[0].connected_components()[1])
            (lowerbound1,Zp1,why1) = msr_bounds(H[1])
            (lowerbound2,Zp2,why2) = msr_bounds(H[2])
            lowerbound = lowerbound1 + lowerbound2 - 1
            why = why1 + why2
            return (lowerbound,why)
		else:
            return 'No cut vertex'
    else:
			return 'No cut vertex'

def dup_vertex_machine(H,just_verts=True):
	"""
	Returns all the duplicate vertices of a graph H or runs
	msr_bounds on the graph without the duplicate vertices.

    INPUT:
        H			--	the graph
        just_verts	--	request to return only the duplicate vertices
						of the graph H (default True)

    OUTPUT:
        If just_verts is True (default), return a list of the sets of 
		duplicate vertices.
		        
       Else, delete all but one copy of each of the "duplicate vertices,"
	   count the number of deletions and return 
	   msr_bounds(graph with deltions, # of deletions).

    EXAMPLES:
		sage: dup_vertex_machine(graphs.CompleteGraph(5))
		[{0, 1, 2, 3, 4}]
		sage: dup_vertex_machine(atlas_graphs[130])
		[{1, 2}, {4, 5}]
		sage: dup_vertex_machine(graphs.PathGraph(5))
		'No duplicate vertices in the graph'
	"""
    
    G = Graph(H.adjacency_matrix())
    adjacencies = []
    for i in G.vertex_iterator():
        adjacent = [i]
        for j in G.vertex_iterator():
            if (min(i,j),max(i,j),None) in G.edges():
                adjacent.append(j)
        adjacent = Set(adjacent)
        adjacencies.append(adjacent)
    
    duplicate_verts = []
    for i in range(0,len(adjacencies)):
        duplicates_of_i = Set([i])
        for j in range(0,len(adjacencies)):
            if j != i:
                A = adjacencies[i]
                B = adjacencies[j]
                if A.intersection(B) == A and A.intersection(B) == B:
                    duplicates_of_i = duplicates_of_i.union(Set([j]))
        if duplicates_of_i.cardinality() > 1:
            duplicate_verts.append(duplicates_of_i)
            
    duplicate_verts = list(Set(duplicate_verts))
    
    if len(duplicate_verts) == 0:
        if just_verts == True:
            return 'No duplicate vertices in the graph'
        else:
            return False
    else: 
        if just_verts == True:
            return duplicate_verts
        else:
            count = 0
            for i in range(0,len(duplicate_verts)):
                for j in range(0,len(duplicate_verts[i])-1):
                    G.delete_vertex(list(duplicate_verts[i])[j])
                    count = count+1
            return (msr_bounds(G),count)

def conn_component_machine(H,just_components=True):
	"""
	Returns connected components of a graph or 
	runs msr_bounds on all connected components.

    INPUT:
        H				--	the graph
        just_components	--	request to return only the connected
							components of the graph H (default True)

    OUTPUT:
		If just_components is True (default), return list of graphs
		which are the connected components of H.
		Else, return sum of msr_bounds run on each component.

    EXAMPLES:
		sage: 
	
	"""
	G = Graph(H.adjacency_matrix())
	if G.is_connected() == True:
		if just_components == True:
			return G
		else:
			return msr_bounds(G)
	else:
		B = G.connected_components()
		components = []
		lowerbound = 0
		Zp = 0
		why = ''
		for i in range(0,G.connected_components_number()):
			Geye=G.copy()
            #
            #
            #Delete vertices from other components
            #
            for k in range(0,G.connected_components_number()):
                if k != i:
                    Geye.delete_vertices(B[k])
					
			if just_components == True:
				components.append(G)
			else:
				(lowerboundeye,Zpeye,whyeye) = msr_bounds(Geye)
				lowerbound = lowerbound + lowerboundeye
				Zp = Zp + Zpeye
				why = why + whyeye
		
		if just_components == True:
			return components
		else:
			return (lowerbound,Zp,why)

def msr_bounds(H):
	"""
	Compute bounds for the minimum semidefinite rank
	of a simple graph H

    INPUT:
        H				--	the graph

    OUTPUT:

    EXAMPLES:
	
	"""

    G = Graph(H.adjacency_matrix())
    
    Zp = Zplus_gen(G)[0]
    ord = G.num_verts()
    lowerbound = 0
    Mp_found = False
    why = ''
	
	#Deal with disconnected graphs
	if G.is_connected() == False:
		(lowerbound,Zp,whydis) = conn_component_machine(G,False)
		why = 'Disconnected: (' + whydis + ')'
    else:
		#If Z+(G) <= 3, M+(G) = Z+(G)
		if Zp <= 3:
			lowerbound = Zp
			Mp_found = True
			why = why + 'Z+(G) <= 3, '
		#M+(G) = 1 <=> Z+(G) = 1 and M+(G) = 2 <=> Z+(G) = 2, so M+(G) >= 3
		else: 
			lowerbound = 3
			case = 'Z+(G) >= 4'
    
		#Z+(G) = |G| - 1 => G is complete, i.e. M+(G) = Z+(G)
		if Zp == ord - 1 and Mp_found == False:
			lowerbound = Zp
			Mp_found = True
			why = why + 'G is a complete graph, '
		
		#If G is chordal, M+(G) = Z+(G)
		if Mp_found == False and G.is_chordal() == True:
			lowerbound = Zp
			Mp_found = True
			why = why + 'G is chordal, '
		
		#Compute k(G)
		if Mp_found == False:
			kappa = int(G.vertex_connectivity())
			#k(G)<=M+(G)<=Z+(G)
			if kappa == Zp:
				lowerbound = kappa
				Mp_found = True
				why = why + 'kappa(G) = Z+(G), '
			else:
				if kappa > lowerbound:
					lowerbound = kappa
					case = 'kappa(G) <= M+(G)'
			
		#Check for duplicate vertices        
		if Mp_found == False:
			check_dup_verts = dup_vertex_machine(G,False)
			if check_dup_verts != False:
				lowerbound = max(check_dup_verts[0][0] + check_dup_verts[1],lowerbound)
				if lowerbound == Zp:
					Mp_found = True
					why = why + 'Dup-vertex: (' + check_dup_verts[0][2] + '), ' 
    
		#Run cut-vertex algorithm    
		if Mp_found == False and kappa == 1:
			(cut_v_lowerbound,cut_v_why) = cut_vertex_machine(G)
			lowerbound = max(cut_v_lowerbound,lowerbound)
			why = why + 'Cut-vertex: (' + cut_v_why + '), '
			if lowerbound == Zp:
				Mp_found = True
    
		#|G|-ccn(G)<=M+(G)<=Z+(G)
		if Mp_found == False:
			lowerbound = max(kappa,lowerbound)
			ccn = len(edge_clique_cover_minimum(G))
			if ccn == ord - Zp:
				lowerbound = ord - Zp
				Mp_found = True
				why = why + 'ccn(G) = |G| - Z+(G), '
			else:
				if G.num_verts() - ccn > lowerbound:
					lowerbound = G.num_verts() - ccn
					case = 'ccn(G) <= |G| - M+(G)'
				
		#Check tree size conditions        
		if Mp_found == False: 
			(works,ts) = msr_is_ts(G)
			if works == 'YES':
				lowerbound = ord - ts + 1
				Mp_found = True
				why = why + 'mr+(G) = ts(G) - 1, '
				
		#Give reasons for lowerbound		
        if Mp_found == False:
			why = why + case + ', '
			    
    return (lowerbound,Zp,why)
			