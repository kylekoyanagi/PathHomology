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
            V = string.(V)
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

    function buildAdjacencyList(V,E)
        adj = [[] for _ in 1:length(V)]
        for e in E
            addEdge!(adj, e[1], e[2])
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
end