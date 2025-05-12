module pathcomplex
#==============================================================================================================
cartesianProduct input: 
(1) set1: a set 
(2) set2: a set

cartesianProcudt Returns: 
(1) the cartesian product of set1Xset2
==============================================================================================================#
    function cartesianProduct(set1,set2)
        result = []
        for i in set1    
            for j in set2

                # try to fix this condition to make it more clean!!!
                if typeof(i) != Int64
                    arr = []
                    for ii in i 
                        push!(arr,ii)
                    end 
                    push!(arr,j)
                    push!(result,arr)
                else
                    push!(result, [i,j])
                end
            end
        end
        return result
    end

#==============================================================================================================
cartesianProductN input: 
(1) set1: a set 
(2) set2: a set
(3) n: an integer, which specifies the number of times you want to take the cartesian product

cartesianProcudtN Returns: 
(1) the cartesian product n times
==============================================================================================================#
    function cartesianProductN(set1,set2,n)
        product = cartesianProduct(set1,set2)

        for i in 1:(n-2) 
            product = cartesianProduct(product,set2)
        end

        return product
    end


#==============================================================================================================
isSubset input: 
(1) set1: a set 
(2) set2: a set

isSubset Returns: 
(1) true, if set1 is contained in set2; false, otherwise
==============================================================================================================#
    function isSubset(set1,set2)
        count = 0
        len = length(set1)
        for j in 1:length(set2)
            count = 0
            for i in 1:length(set1)
                if length(set1) < (j + count)
                    break
                end

                if set1[i] != set2[j+count]
                    break
                else
                    count = count + 1
                    if count == len
                        return true
                    end
                end
            end
        end
        return false
    end


    function isSubset2(set1, set2)
        if length(set1) <= length(set2)
            for j in 1:length(set2)
                position = 0
                for i in 1:length(set1)

                    if length(set2) < (j + position)
                        break
                    end

                    if set1[i] != set2[j + position]

                        break
                    else
                        position = position + 1
                    end  
                end
                if position == length(set1)
                    return true
                end
            end
            return false
        else 
            return false
        end
    end

#==============================================================================================================
maximalEdgeSet input: 
(1) E: the edge set of a hypergraph

maximalEdgeSet Returns: 
(1) the maximal edge set
==============================================================================================================#
    function maximalEdgeSet(E)
        maximalE = []

        for i in 1:length(E)
            e = E[i]
            maximal = true
            for j in 1:length(E)
                if i == j 
                    continue
                else
                    testE = E[j]
                    if isSubset2(e, testE) == true
                        maximal = false
                        break
                    end
                end
            end
            if maximal == true
                push!(maximalE, e)
            end
        end
        return maximalE
    end


    function buildpathcomplex(H,q,n)
        P = []
        vertices = H[1]
        edges = H[2]
        push!(P,vertices)

        for i in 2:n
            Pi = []
            Vi = cartesianProductN(vertices,vertices,i)
            for p in Vi
                validPath = true
                for j in 1:length(p)
                    if  length(p) < (j+q-1)
                        break
                    end 
                    slice = p[j:(j+q-1)]
                    subset = false
                    for e in edges
                        if isSubset2(slice,e) == true
                            subset = true
                        end
                    end
                    if subset == false
                        validPath = false
                        break
                    end
                end
                if validPath == true
                    push!(Pi,p)
                end
            end
            push!(P,Pi)
        end
        return P
    end

#==============================================================================================================
removeSelfLoops input: 
(1) P: a set of paths of length 2

removeSelfLoops Returns: 
(1) p: the set of path of length 2 such that i_0 != i_1
==============================================================================================================#
    function removeSelfLoops(P)
        p = []
        for i in 1:length(P)
            if P[i][1] != P[i][2]
                push!(p, P[i])
            end
        end
        return p
    end


    function removeNselfLoops(P)
        p = []
        for i in 1:length(P)
            regular = true
            path = P[i]
            for j in 1:(length(path)-1)
                if path[j] == path[j+1]
                    regular = false
                    break
                end
            end
            if regular == true
                push!(p, path)
            end
        end
        return p
    end 
#==============================================================================================================
buildpathcomplexV2 input: 
(1) H: hypergraph 
(2) q: integer, 
(3) n: ineger, up to the length of paths to take

buildpathcomplexV2 Returns: 
(1) P: path complex, specifically P^q(H).
==============================================================================================================#
    function buildpathcomplexV2(H,q,n)
            P = []
            vertices = H[1]
            edges = H[2]
            maximalEdges = maximalEdgeSet(edges)
            push!(P,string.(vertices))

            for i in 2:n
                Pi = []
                Vi = cartesianProductN(vertices,vertices,i)
                for p in Vi
                    validPath = true
                    if length(p) < q
                        subset = false
                        for e in maximalEdges
                            if isSubset2(p,e) == true
                                subset = true
                            end
                        end
                        if subset == false
                            validPath = false
                        end

                    else
                        for j in 1:length(p)
                            if  length(p) < (j+q-1)
                                break
                            end 
                            slice = p[j:(j+q-1)]
                            subset = false
                            for e in maximalEdges
                                if isSubset2(slice,e) == true
                                    subset = true
                                end
                            end
                            if subset == false
                                validPath = false
                                break
                            end
                        end
                    end
                    if validPath == true
                        push!(Pi,string.(p))
                    end
                end

               #= if i == 2
                    Pi = removeSelfLoops(Pi)
                end=#
                Pi = removeNselfLoops(Pi)
                push!(P,Pi)
            end
            return P
        end
    end

