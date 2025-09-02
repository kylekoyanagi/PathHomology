module graph_preprocess
    
    function relabel(V,E) 
        newV = []
        newE = [] 
        vertD = Dict()
        for i in 1:length(V)
            push!(newV, i)
            vertD[V[i]] = i
        end
        for e in E

            push!(newE, [vertD[e[1]],vertD[e[2]]])
        end
        return newV, newE
    end   

    function relabelGraph(V,E)
        if !(typeof(V) == Vector{Int64})
            V = V
            V, E = relabel(V,E)
        else
            if minimum(V) < 1
                V = string.(V)
                V, E = relabel(V,E)
            end
        end
        return V,E
    end

    function addEdge!(adj,u,v)
        push!(adj[u], v)
    end

    function addWeightedEdge!(adj,u,v,w)
        push!(adj[u], (v,w))
    end

    function buildAdjacencyList(V,E)
        adj = [[] for _ in 1:length(V)]
        for e in E
            addEdge!(adj, e[1], (e[2]))
        end
        return adj
    end

    function buildWeightedAdjacencyList(V,E)
        adj = [[] for _ in 1:length(V)]
        for e in E
            addWeightedEdge!(adj, e[2][1], e[2][2], e[1])
        end
        return adj
    end

    function removeHangingVertex(V,E, adj)
        for i in 1:length(adj)
            aList = adj[i]
            if length(aList) == 1
                contained = false
                for list in adj
                    if i in list
                        contained = true 
                        break
                    end
                end
                if contained == false
                    deleteat!(V, findall(x->x==i,V))
                    deleteat!(E, findall(x->x==[i,aList[1]],E))
                    return V,E
                end
            end

            if length(aList) == 0
                count = 0
                index = -1
                for j in 1:length(adj)
                    list = adj[j]
                    if i in list 
                        index = j
                        count += 1
                    end
                end
                if count == 1
                    deleteat!(V, findall(x->x==i,V))
                    deleteat!(E, findall(x->x==[index,i],E))
                    return V,E
                end
            end
        end
        return V,E
    end

    function removeDegreeNVertex(V,E, adj)
        for i in 1:length(adj)

            if isempty(adj[i]) == true
                list = []
                for j in 1:length(adj)
                    if !(i == j)
                        if i in adj[j]
                            push!(list, j)
                        end
                    end
                end

                for l in list 
                    count = 0 
                    for k in list 
                        if !(i == k)
                            if l in adj[k]
                                count += 1
                            end
                        end
                    end
                    if count == (length(list)-1)
                        deleteat!(V, findall(x->x==i,V))
                        for bi in list
                            deleteat!(E, findall(x->x==[bi,i],E))
                        end
                        break
                    end
                end

            else
                contained = false
                for alist in adj
                    if i in alist 
                        contained = true
                        break
                    end
                end

                if contained == false
                    iList = adj[i]
                    for v in iList 
                        count = 0
                        vList = adj[v]
                        for j in iList 
                            if j in vList 
                                count += 1
                            end
                        end
                        if count  == (length(iList) - 1)
                            deleteat!(V, findall(x->x==i,V))
                            for bi in iList
                                deleteat!(E, findall(x->x==[i,bi],E))
                            end
                            break
                        end
                    end
                end
            end
        end
        return V,E
    end

    function isSemiEdge(a,b,c,adj)
        for d in adj[c]
            if !(a == d)
                if b in adj[d]
                    return true
                end
            end
        end
        return false
    end

    function removeDegree1plus1Vertex(V,E,adj)
        for a in V 
            if length(adj[a]) == 1
                b = adj[a][1]
                count = 0
                index = -1
                for c in V 
                    if a in adj[c]
                        count += 1
                        index = c 
                    end
                end
                if count == 1
                    if b in adj[index]
                        deleteat!(V, findall(x->x==a,V))
                        deleteat!(E, findall(x->x==[a,b],E))
                        deleteat!(E, findall(x->x==[index,a],E))
                        break
                    else
                        if isSemiEdge(a,b,index,adj) == true
                            deleteat!(V, findall(x->x==a,V))
                            deleteat!(E, findall(x->x==[a,b],E))
                            deleteat!(E, findall(x->x==[index,a],E))   
                            break 
                        end
                    end
                end
            end
        end
        return V,E
    end

    function cleanGraph(V,E)
        V, E = relabelGraph(V,E)
        adj = buildAdjacencyList(V,E)
        remove = true 
        n = length(V)
        hanging = false
        degree_n = false
        deg_1plus1 = false

        for i in 1:length(V)
            V,E = removeHangingVertex(V,E,adj)
            if length(V) < n
                hanging = true
                V,E = relabel(V,E)
                adj = buildAdjacencyList(V,E)
                n = length(V)

            else 
                hanging = false
            end

            V,E = removeDegreeNVertex(V,E,adj)
            if length(V) < n
                degree_n = true
                V,E = relabel(V,E)
                adj = buildAdjacencyList(V,E)
                n = length(V)
            else 
                degree_n = false
            end
            V,E = removeDegree1plus1Vertex(V,E,adj)
            if length(V) < n
                deg_1plus1 = true
                V,E = relabel(V,E)
                adj = buildAdjacencyList(V,E)
                n = length(V)
            else 
                deg_1plus1 = false
            end

            if (degree_n == false) & (hanging = false) & (deg_1plus1 == false)
                break
            end
        end
        return V,E
    end

    function makeSymmetricReflexive(A)
        B=copy(A)
        for a in A
            if !( [a[2],a[1]] in B)
                push!(B,[a[2],a[1]] )
            end
        end
        return B
    end
    
    function cartesian_product(arrays) 
        # Base case: if arrays is empty, return an array with an empty tuple
        if isempty(arrays)
            return [()]
        end
        
        # Get the first array and the rest of the arrays
        first_array, rest_arrays = arrays[1], arrays[2:end]
        
        # Recursive call to get the Cartesian product of the rest of the arrays
        rest_product = cartesian_product(rest_arrays)
        
        # Combine each element of the first array with each tuple from the rest product
        result = [(x, rest...) for x in first_array for rest in rest_product]
        
        return result
    end
    
    function multBoxProd(graphs) #input list of graphs
    
        vert=cartesian_product([G.vertices for G in graphs])
        edge=[]
        
        for v in vert
            
            for w in vert
                tot = 0 # check if its an edge in the product
                
                for i = 1:length(v)
                    
                    if !(v[i] == w[i]) && ( [v[i],w[i]] in graphs[i].edges) #check if theye connected in g_i
                        tot += 1
                    elseif !([v[i],w[i]] in graphs[i].edges)
                        tot += 2
                    end
                    
                end
                if tot <2
                    push!(edge,[v,w])
                end
            end
        end
    
        return(digraph(vert,edge))
    
    end
    
    function cartesian_product(arrays) 
        # Base case: if arrays is empty, return an array with an empty tuple
        if isempty(arrays)
            return [()]
        end
        
        # Get the first array and the rest of the arrays
        first_array, rest_arrays = arrays[1], arrays[2:end]
        
        # Recursive call to get the Cartesian product of the rest of the arrays
        rest_product = cartesian_product(rest_arrays)
        
        # Combine each element of the first array with each tuple from the rest product
        result = [(x, rest...) for x in first_array for rest in rest_product]
        
        return result
    end
    
    #Allows you to relate two vertices in a graph
    function quotient_vertices(G,v1,v2)
        new_vert=[]
        new_edges=[]
    
        for v in G.vertices
            if !(v==v2)
                push!(new_vert,v)
            end
        end
    
        for i in G.edges
            if !(i[1]==v2 || i[2]==v2)
                push!(new_edges,i)
            end
            if i[1]==v2 && !(i[2]==v2) && !(i[2]==v1)
                push!(new_edges,[v1,i[2]])
            end
            if !(i[1]==v2) && i[2]==v2 && !(i[1]==v1)
                push!(new_edges,[i[1],v1])
            end
        end
        return digraph(unique!(new_vert),unique!(new_edges))
    end
    
    #Allows you to relate two vertices in a graph
    function quotient_vertices(G,v1,v2)
        new_vert=[]
        new_edges=[]
    
        for v in G.vertices
            if !(v==v2)
                push!(new_vert,v)
            end
        end
    
        for i in G.edges
            if !(i[1]==v2 || i[2]==v2)
                push!(new_edges,i)
            end
            if i[1]==v2 && !(i[2]==v2) && !(i[2]==v1)
                push!(new_edges,[v1,i[2]])
            end
            if !(i[1]==v2) && i[2]==v2 && !(i[1]==v1)
                push!(new_edges,[i[1],v1])
            end
        end
        return digraph(unique!(new_vert),unique!(new_edges))
    end
end