module persistence
    include("homology.jl")
    include("graph_preprocess.jl")
    include("smith.jl")
    using SparseArrays;
    using LinearAlgebra;
    using .SNF
    using .homology
#=
NEED TO DO: 
    If not taking the entire filtration, 
    need to consider paths that are not allowed after taking the differential, 
    if a path is not allowed, then add at infinity 



=#
    function filtration(G,D,F)
        edges = G[2]
        filteredGraphs = []
        for f in F
            filteredV = []
            filteredE = []
            for e in edges
                vertexIn = e[1]
                vertexOut = e[2]
                d = D[vertexIn,vertexOut]
                if d <= f
                    push!(filteredE, e)
                    push!(filteredV, vertexIn)
                    push!(filteredV, vertexOut)
                end
            end
            if isempty(filteredV) == false 
                push!(filteredGraphs,[filteredV,filteredE])
            end
        end
        return filteredGraphs
    end

    function dijkstrasAlg(V,E,s)
        adj = graph_preprocess.buildWeightedAdjacencyList(V,E)
        dist = [Inf for i in 1:length(V)] 
        visited = [false for i in 1:length(V)]
        path = [0 for i in 1:length(V)]

        for i in 1:length(V)
            path[i] = -1
        end

        dist[s] = 0
        path[s] = -1
        current = s

        sett = Set()

        #for k in 1:(length(V))
        while true
            visited[current] = true
            for i in 1:length(adj[current])
                v = adj[current][i][1]
                if visited[v]
                    continue 
                end 

                push!(sett, v)
                alt  = dist[current] + adj[current][i][2]

                if alt < dist[v]
                    dist[v] = alt
                    path[v] = current 
                end
            end

            if current in sett 
                delete!(sett,current)
            end

            if length(sett) == 0
                break
            end

            minDist = Inf
            index = 0 

            for a in sett 
                if dist[a] < minDist
                    minDist = dist[a]
                    index = a 
                end
            end
            current = index
        end
        return dist
    end

    function shortestPathFiltration(G,F)
        wghtE = G[2]
        V = G[1]
        filteredE = []
        for v in V 
            dist = dijkstrasAlg(V,wghtE, v)
            for i in 1:length(dist)
                if !(i == v) & !(dist[i] == Inf)
                    push!(filteredE, [dist[i],[Int64.(v),Int64.(i)]])
                end
            end
        end
        sortedFilteredE = sort(filteredE)
        sortedE = []
        sortedW = []
        for i in 1:length(sortedFilteredE)
            push!(sortedE,sortedFilteredE[i][2])
            push!(sortedW, sortedFilteredE[i][1])
        end
        F = sort(F)
        filteredGraph = []
        for filtration in F 
            indices = findall(x->x<=filtration,sortedW)
            push!(filteredGraph, [V, [sortedW[indices], sortedE[indices]]])
        end
        return filteredGraph
    end

    # regularPaths builds all regular paths up to paths of length n
    function regularPaths(V,n)
        RP = []
        push!(RP, V)
        for i in 1:n
            RP_1 = RP[i] 
            paths = [] 
            for rp in RP_1 
                for v1 in V
                    if i == 1  
                        if !(v1 == rp)
                            push!(paths, (Inf,[rp, v1]))
                        end
                    else
                        temp = []
                        if !(v1 == rp[2][length(rp)])
                            temp = copy(rp[2]) 
                            push!(temp, v1)
                            push!(paths, (Inf,temp))
                        end
                    end
                end
            end
            push!(RP, paths)
        end
        return RP
    end

    # regularPaths builds all regular paths up to paths of length n
    function regularPaths2(V,n)
        RP = []
        push!(RP, V)
        for i in 1:n
            RP_1 = RP[i] 
            RP_1Keys = keys(RP_1)
            paths = Dict() 
            for rp in RP_1Keys 

                for v1 in V
                    if i == 1  
                        if !(v1 == rp)
                            paths[[rp,v1]] = Inf
                        end
                    else
                        temp = []
                        if !(v1 == rp[length(rp)])
                            temp = copy(rp) 
                            push!(temp, v1)
                            paths[temp] = Inf
                        end
                    end
                end
            end
            push!(RP, paths)
        end
        return RP
    end

    # computeAllowTime computes the allow time of all allowed paths of a filtered graph 
    # it suffices to use the weights from the last filtration (the biggest filtration value) as all smaller are contained in this graph
    function computeAllowTime(G,F,n)
        filteredGraphs = shortestPathFiltration(G,F)
        vertices = filteredGraphs[1][1]
        A = [] 
        numFiltrations = length(F)
        push!(A, vertices)
        for i in 1:n 
            paths = Dict()
            if i == 1 
                filteredGraph = filteredGraphs[numFiltrations][2] 
                weights = filteredGraph[1]
                edges = filteredGraph[2]
                for j in 1:length(weights)
                    paths[edges[j]] = weights[j]
                end
            else
                edgeKeys = keys(A[2])
                for edge in edgeKeys
                    pathKeys = keys(A[i])
                    for path in pathKeys
                        if edge[1] == path[length(path)] # if the start of an edge equals the end of the path
                            ap = vcat(path,[edge[2]])
                            paths[ap] = max(A[2][edge], A[i][path])
                        end
                    end
                end
            end
            push!(A, paths)
        end
        return A
    end

    # updateRegularPaths updates the weights of regular paths if it is an allowed path 
    function updateRegularPaths(A,RP,n)

        for i in 2:n 
            Ai = A[i]
            AiKeys = keys(Ai)
            RPi = RP[i]

            for aPath in AiKeys 
                Aweight = Ai[aPath]
                RPi[aPath] = Aweight                 
            end
        end
        return RP
    end

    # dictToList turns a dicitonary into a list of the form [[key1, value1], ..., [keyn, valuen]]
    function dictToList(D)
        list = [] 
        dKeys = collect(keys(D))

        for k in dKeys 
            push!(list, [D[k], k])
        end
        return list
    end

    # Extended Euclidean Algorithm 
    function gcdExt(a, b)
        if a == 0
            return b, 0, 1
        end 
        m = mod(b,a)
        gcd, x1, y1 = gcdExt(m, a)
    
        x = y1 - floor(b/a)*x1
        y = x1
        return gcd, x, y
    end

    function getRange(wghts)
        ub = []
        for i in 1:(length(wghts)-1)
            if !(wghts[i] == wghts[i+1])
                push!(ub, i)
            end
        end
        range = []

        if isempty(ub) == false
            for j in 1:length(ub)
                if j == 1 

                    if length(ub) == 1
                        push!(range, [1,ub[1]])
                        push!(range, [ub[1]+1, length(wghts)])
                    else
                        push!(range, [1,ub[1]])
                    end

                elseif j == length(ub) 
                    if ub[j] == length(wghts)
                        push!(range, [ub[j-1]+1, ub[j]])
                    else
                        push!(range, [ub[j-1]+1, ub[j]])
                        push!(range, [ub[j]+1, length(wghts)])
                    end

                else
                    push!(range,[ub[j-1]+1, ub[j]])
                end
            end
        else
            push!(range, [1, length(wghts)])
        end
        return range
    end

    # sort regular paths lexographically within same path allow time
    function pathSort(sortedRP)
        wghts = [p[1] for p in sortedRP]
        wghtRng = getRange(wghts)

        sortedPaths = []
        for rng in wghtRng 
            s = sort!(sortedRP[rng[1]:rng[2]], by = x -> x[2])
            sortedPaths = vcat(sortedPaths, s)
        end

        return sortedPaths 
    end

    #buildMatrix builds a matrix with RP as the domain and RP_1 as the codomain 
    function buildMatrix(RP,RP_1,P_1)
        sortedRP = sort(dictToList(RP))
        #println("SORTEDRP; ", sortedRP)
        wghts = [p[1] for p in sortedRP]

        weightRange = getRange(wghts)
        ET = Dict()
        M = []
        if P_1 == 1
            
            for i in 1:length(sortedRP)
                rpInfo = sortedRP[i]
                rp = rpInfo[2]
                boundaryRP = homology.dV2(rp,2,1)
                col = zeros(length(RP_1))
                boundaryKeys = keys(boundaryRP)

                for boundaryPath in boundaryKeys 
                    index = homology.getIndex(RP_1,boundaryPath[1])

                    col[index] = col[index] + boundaryRP[boundaryPath]
                end
                push!(M, col)
            end
 
        else
            #sortedRP_1 = sort(RP_1, rev = true)
            sortedRP_1 = sort(dictToList(RP_1), rev=true)
            #sortedRP_1 = pathSort(sortedRP_1)
            println("sortedRP_1: ", sortedRP_1)

            pathList = [] 

            for i in 1:length(sortedRP_1)
                info = sortedRP_1[i] 
                p = info[2]
                push!(pathList, p)
            end

            for i in 1:length(sortedRP)
                rpInfo = sortedRP[i]
                rp = rpInfo[2]
                boundaryRP = homology.dV2(rp,P_1+1,1)
                col = zeros(length(pathList))
                boundaryKeys = keys(boundaryRP)
                maxindex = Inf
                for boundaryPath in boundaryKeys 
                    index = homology.getIndex(pathList,boundaryPath)
                     maxindex = Int.(min(maxindex, index))
                    if index == -1 
                        #println("index is -1: ", (boundaryPath, rp))
                    end
                    col[index] = col[index] + boundaryRP[boundaryPath]
                end
                #ET[sortedRP_1[maxindex][2]] = sortedRP_1[maxindex][1]
                ET[rp] = sortedRP_1[maxindex][1]
                push!(M, col)
            end
        end
        return sparse(reduce(hcat,M)), weightRange, ET
    end

    # Swap Columns
    function swapcol!(x,i,j)
        for k in axes(x,1) 
        idata = x[k,i]
        x[k,i] = x[k,j]
        x[k,j] = idata
        end
    end

    # columnReduceMatrix performs left to right column Gaussian Elimination to reduce the matrix 
    # a transformation matrix T is calculated to get the linear combinations of the omega basis elements 
    function columnReduceMatrix(M, wghtRange)
        #= 
            COLUMN REDUCTION STEPS: 
                (1) Find Pivot in each column, the first nonzero row
                (2) Clear Column

        =#
        m,n = size(M)
        T = sparse(1I, n, n)

        for i in 1:n
            rangeIndex = 0
            for j in 1:length(wghtRange)
                if (i <= wghtRange[j][2]) & (i >= wghtRange[j][1])
                    rangeIndex = j 
                    break 
                end
            end

            possiblePivots = []

            for k in i:wghtRange[rangeIndex][2]
                possibleCol = M[:,k]
                if isempty(findnz(possibleCol)[1]) == false
                    possiblePivot = minimum(findnz(possibleCol)[1])
                    push!(possiblePivots, (possiblePivot, k))
                end
            end

            if isempty(possiblePivots) == false
                smallestPivot = sort(possiblePivots)[1][2]

                if !(smallestPivot == i)
                    swapcol!(M,i, smallestPivot)
                    swapcol!(T,i, smallestPivot)
    
                end
            end

            col = M[:,i] 
            if isempty(findnz(col)[1]) == false
                pivot = minimum(findnz(col)[1])
                nonzeroCols = findnz(M[pivot, : ])[1]

                for index in nonzeroCols 
                    if index > i
                        pivotval = M[pivot, i]
                        nonzeroEntry = M[pivot, index]
                        if pivotval == 1
                            M[:, index] = M[:, index] - nonzeroEntry*M[:, i]
                            T[:, index] = T[:, index] - nonzeroEntry*T[:, i]

                        elseif (mod(nonzeroEntry,pivotval) == 0)
                            factor = nonzeroEntry / pivotval
                            
                            M[:, index] = M[:, index] - factor*M[:, i]
                            T[:, index] = T[:, index] - factor*T[:, i]

                        else
                            d, x, y = gcdExt(pivotval, nonzeroEntry)

                            # scale col
                            M[:, index] = y*M[:, index]
                            T[:, index] = y*T[:, index]

                            # col addition to get M[pivot, index] = d
                            M[:, index] = M[:, index] + x*M[:, i]
                            T[:, index] = T[:, index] + x*T[:, i]

                            factor = pivotval / d
                            # scale col d*factor
                            M[:, index] = factor*M[:, index]
                            T[:, index] = factor*T[:, index]

                            # col subtraction
                            M[:, index] = M[:,index] - M[:, i]  
                            T[:, index] = T[:,index] - T[:, i]        
      
                        end
                    end
                end
                dropzeros!(M)
            end
            dropzeros!(M)
        end
        dropzeros!(T)
        return M,T
    end

    function findPivot(M,t)
        tRow = M[t,:]
        nonzeroCols = findnz(tRow)[1]

        for col in nonzeroCols 
            if col >= t
                return col
            end
        end
        return -1
    end

    function findPivot2(M,t, i)
        tRow = M[t,:]
        nonzeroCols = findnz(tRow)[1]

        for col in nonzeroCols 
            if col >= i
                return col
            end
        end
        return -1
    end

    function findnzcols(M)
        ncols = size(M)[2]
        nzcols = []
        pivotcols = []
        for i in 1:ncols 
            col = M[:,i]
            colEntries = findnz(col)[1]
            if isempty(colEntries) == true
                push!(nzcols, i)
            else
                push!(pivotcols, i)
            end
        end
        return nzcols, pivotcols
    end

        # columnReduceMatrix performs left to right column Gaussian Elimination to reduce the matrix 
    # a transformation matrix T is calculated to get the linear combinations of the omega basis elements 
    function columnReduceMatrix2(M)
        #println(M)
        #= 
            COLUMN REDUCTION STEPS: 
                (1) Find Pivot in each column, the first nonzero row
                (2) Clear Column

        =#
        m,n = size(M)
        column_tracker = collect(1:n)
        T = sparse(1I, n, n)
        for i in 1:n
            if i > m
                break 
            end
 
            index = i 
            pivot = findPivot2(M, i, i)
            cont = true
            if pivot == -1
                while pivot == -1
                    if index > m
                        #println((index,m)) 
                        cont = false
                        break
                    end
                    pivot = findPivot2(M, index, i)
                    index = index + 1
                end
            end
            # if cont. is false then the rest of the matrix is all zeros 
            if cont == false
                continue
            end
            column_tracker[[pivot,i]] = column_tracker[[i,pivot]]
            swapcol!(M,pivot, i)
            swapcol!(T, pivot, i)
            col = M[:,i] 
            if isempty(findnz(col)[1]) == false
                pivot = minimum(findnz(col)[1])
                nonzeroCols = findnz(M[pivot, : ])[1]

                for index in nonzeroCols 
                    if index > i
                        pivotval = M[pivot, i]
                        nonzeroEntry = M[pivot, index]
                        if pivotval == 1
                            M[:, index] = M[:, index] - nonzeroEntry*M[:, i]
                            T[:, index] = T[:, index] - nonzeroEntry*T[:, i]

                        elseif (mod(nonzeroEntry,pivotval) == 0)
                            factor = nonzeroEntry / pivotval
                            
                            M[:, index] = M[:, index] - factor*M[:, i]
                            T[:, index] = T[:, index] - factor*T[:, i]

                        else
                            d, x, y = gcdExt(pivotval, nonzeroEntry)

                            # scale col
                            M[:, index] = y*M[:, index]
                            T[:, index] = y*T[:, index]

                            # col addition to get M[pivot, index] = d
                            M[:, index] = M[:, index] + x*M[:, i]
                            T[:, index] = T[:, index] + x*T[:, i]

                            factor = pivotval / d
                            # scale col d*factor
                            M[:, index] = factor*M[:, index]
                            T[:, index] = factor*T[:, index]

                            # col subtraction
                            M[:, index] = M[:,index] - M[:, i]  
                            T[:, index] = T[:,index] - T[:, i]        
    
                        end
                    end
                end
                dropzeros!(M)
            end
            dropzeros!(M)
        end

        dropzeros!(T)
        dropzeros!(M)
        #println("COL TRACKER: ", column_tracker)
        nzcols, pivotcols = findnzcols(M)
        return M,T, nzcols, pivotcols, column_tracker
    end

    function columnReduceMatrix3(M)
        #println(M)
        #= 
            COLUMN REDUCTION STEPS: 
                (1) Find Pivot in each column, the first nonzero row
                (2) Clear Column

        =#
        m,n = size(M)
        column_tracker = collect(1:n)
        T = sparse(1I, n, n)
        for i in 1:n
            index = i
            col = M[:,i] 
            if isempty(findnz(col)[1]) == false
                pivot = minimum(findnz(col)[1])
                nonzeroCols = findnz(M[pivot, : ])[1]

                for index in nonzeroCols 
                    if index > i
                        pivotval = M[pivot, i]
                        nonzeroEntry = M[pivot, index]
                        if pivotval == 1
                            M[:, index] = M[:, index] - nonzeroEntry*M[:, i]
                            T[:, index] = T[:, index] - nonzeroEntry*T[:, i]

                        elseif (mod(nonzeroEntry,pivotval) == 0)
                            factor = nonzeroEntry / pivotval
                            
                            M[:, index] = M[:, index] - factor*M[:, i]
                            T[:, index] = T[:, index] - factor*T[:, i]

                        else
                            d, x, y = gcdExt(pivotval, nonzeroEntry)

                            # scale col
                            M[:, index] = y*M[:, index]
                            T[:, index] = y*T[:, index]

                            # col addition to get M[pivot, index] = d
                            M[:, index] = M[:, index] + x*M[:, i]
                            T[:, index] = T[:, index] + x*T[:, i]

                            factor = pivotval / d
                            # scale col d*factor
                            M[:, index] = factor*M[:, index]
                            T[:, index] = factor*T[:, index]

                            # col subtraction
                            M[:, index] = M[:,index] - M[:, i]  
                            T[:, index] = T[:,index] - T[:, i]        
    
                        end
                    end
                end
                #println("Matrix: ", (pivot, M))
                if isempty(M[:,1]) == true
                    println("EMPTYYYYY: ", pivot)
                end
                dropzeros!(M)
            end
            dropzeros!(M)
        end

        dropzeros!(T)
        dropzeros!(M)

        #println("COL TRACKER: ", column_tracker)
        nzcols, pivotcols = findnzcols(M)

        pivotRel = Dict()
        for pivotcol in pivotcols 
            highestentry = minimum(findnz(M[:,pivotcol])[1])
            pivotRel[pivotcol] = highestentry
        end
        return M,T, nzcols, pivotcols, column_tracker, pivotRel
    end

    function columnReduceMatrix4(M)
        #println(M)
        #= 
            COLUMN REDUCTION STEPS: 
                (1) Find Pivot in each column, the first nonzero row
                (2) Clear Column

        =#
        m,n = size(M)
        pivotTracker = zeros(m)
        column_tracker = collect(1:n)
        T = sparse(1I, n, n)
        for i in 1:n
            index = i
            col = M[:,i] 
            if isempty(findnz(col)[1]) == false
                while isempty(findnz(M[:,i])[1]) == false
                    col = M[:,i] 
                    pivot = minimum(findnz(col)[1])
                    if pivotTracker[pivot] == 0
                        pivotTracker[pivot] = i   
                        println(pivotTracker)
                        break 
                    else
                        columnIndex = Int64.(pivotTracker[pivot])
                        
                        pivotval = M[pivot, columnIndex]
                        nonzeroEntry = M[pivot, i]
                        if pivotval == 1
                            M[:, i] = M[:, i] - nonzeroEntry*M[:, columnIndex]
                            T[:, i] = T[:, i] - nonzeroEntry*T[:, columnIndex]

                        elseif (mod(nonzeroEntry,pivotval) == 0)
                            factor = nonzeroEntry / pivotval
                            
                            M[:, i] = M[:, i] - factor*M[:, columnIndex]
                            T[:, i] = T[:, i] - factor*T[:, columnIndex]

                        else
                            d, x, y = gcdExt(pivotval, nonzeroEntry)

                            # scale col
                            M[:, i] = y*M[:, i]
                            T[:, i] = y*T[:, i]

                            # col addition to get M[pivot, index] = d
                            M[:, i] = M[:, i] + x*M[:, columnIndex]
                            T[:, i] = T[:, i] + x*T[:, columnIndex]

                            factor = pivotval / d
                            # scale col d*factor
                            M[:, i] = factor*M[:, i]
                            T[:, i] = factor*T[:, i]

                            # col subtraction
                            M[:, i] = M[:,i] - M[:, columnIndex]  
                            T[:, i] = T[:,i] - T[:, columnIndex]        
                        end
                    end
                    dropzeros!(M)
                end
                #println("Matrix: ", (pivot, M))
                if isempty(M[:,1]) == true
                    println("EMPTYYYYY: ", pivot)
                end
                dropzeros!(M)
            end
            dropzeros!(M)
        end

        dropzeros!(T)
        dropzeros!(M)
        #println("COL TRACKER: ", column_tracker)
        nzcols, pivotcols = findnzcols(M)
        return M,T, nzcols, pivotcols, column_tracker
    end
 
    function getOmegaElts(T, nzcols, RP)
        basisElts = []

        for index in nzcols 
            col = T[:, index] 
            nzrows = findnz(col)[1] 
            lowestEntry = maximum(nzrows)
            basisElt = RP[lowestEntry]
            push!(basisElts, basisElt)
        end
        return basisElts
    end

    function getOmegaElts2(T, nzcols, column_tracker, RP)
        basisElts = []
        for index in nzcols 
            eltIndex = column_tracker[index]
            basisElt = RP[eltIndex]
            push!(basisElts, basisElt)
        end
        return basisElts
    end

    function getPivotElts(M, pivotcols, RP_1)
        pivotElts = []
        for index in pivotcols 
            col = M[:, index] 
            nzrows = findnz(col)[1] 
            highestEntry = minimum(nzrows)
            pivotElt = RP_1[highestEntry]
            push!(pivotElts, pivotElt)
        end
        return pivotElts
    end

    function getPivotElts2(M, pivotcols, RP_1,rev_tracked)
        #println("getPIVOT M: ",M)
        omegaRP_1 = RP_1[rev_tracked]
        pivotElts = []
        RP_1 = sort(RP_1, rev=true)
        #println("sortedOMEGAELTS: ",RP_1)
        rev_tracked = reverse(rev_tracked)
        #println("rev_tracked: ",rev_tracked)
        for index in pivotcols 
            col = M[:, index] 
            #println(col)
            nzrows = findnz(col)[1] 
            highestEntry = minimum(nzrows)
            #println("highest entry: ", highestEntry)

            trackedHighestEntry = rev_tracked[highestEntry]
            pivotElt = RP_1[trackedHighestEntry]
            #println("pivotELT: ", pivotElt)
            #push!(pivotElts, omegaRP_1[highestEntry])

            push!(pivotElts, pivotElt)
        end
        return pivotElts
    end

    # getOmegaBasis from the transformation matrix T, the omega basis is calulated of the form [(weight, [path1, ..., pathn], [sign1, ... signn]), ...]
    function getOmegaBasis(T, RP)  
        sortedRP = sort(dictToList(RP))
        pathList = []
        weightList = []

        for i in 1:length(sortedRP)
            info = sortedRP[i] 
            p = info[2]
            w = info[1]
            push!(pathList, p)
            push!(weightList, w)
        end
        
        newBasis = []
        m,n = size(T) 

        for i in 1:n 
            weight = weightList[i]
            basisElt = []
            signList = []
            col = T[:, i]
            nonzeroRowsInfo = findnz(col)
            nonzeroRows = nonzeroRowsInfo[1]
            nonzeroRowsVals = nonzeroRowsInfo[2]
            for j in 1:length(nonzeroRows)
                row = nonzeroRows[j]
                value = nonzeroRowsVals[j]
                push!(basisElt, pathList[row])
                push!(signList, value)
            end
            push!(newBasis, (weight, basisElt, signList))
        end
        return newBasis
    end

    function buildOmegaBases(regularPaths,n)
        transformationMatrices = []
        reducedMatrices = []
        for i in 2:n
            RP = regularPaths[i]
            RP_1 = regularPaths[i-1]

            M = buildMatrix(RP, RP_1, i-1)

            if ((i-1) == 1)
                Mreduced, T = columnReduceMatrix(M)
                push!(transformationMatrices, T)
                push!(reducedMatrices, Mreduced)
            else    
                transformationMatrix = transformationMatrices[i-2]

                row,col = size(transformationMatrix)
                swapNum = Int64.(floor(col/2)) 
                for j in 1:swapNum 
                    swapcol!(transformationMatrix, j, col+1-j)
                end

                MtimesT = transformationMatrix*M

                Mreduced, T = columnReduceMatrix(MtimesT)

                push!(transformationMatrices, T)
                push!(reducedMatrices, Mreduced)
            end
        end
        return reducedMatrices
    end

    function buildOmegaBases2(regularPaths,n)
        transformationMatrices = []
        reducedMatrices = []
        for i in 2:n
            RP = regularPaths[i]
            RP_1 = regularPaths[i-1]

            M,wghtRng = buildMatrix(RP, RP_1, i-1)
            Mreduced, T = columnReduceMatrix(M,wghtRng)
            push!(transformationMatrices, T)
            push!(reducedMatrices, Mreduced)
        end
        return reducedMatrices
    end

    function calculateEntryTime(pivotElts, pivotdiffElts, basisElts, n)
        #println("PIVOT ELTS: ", pivotElts)
        #println(length(pivotElts))
        entryTime = Dict()
        pivotColRow = Dict()
        if n == 1
            for i in 1:length(pivotElts) 
            
                pivotElt = pivotElts[i]
                diffPivotElt = pivotdiffElts[i]
                pivotColRow[pivotElt] = diffPivotElt
                dweight, dpath = diffPivotElt
                weight, path = pivotElt 
                entryTime[path] = max(weight, 0)
            end
            #println("PCR: ",pivotColRow)
            for i in 1:length(basisElts)
                basisElt = basisElts[i]
                weight, path = basisElt 
                entryTime[path] = max(weight, 0)
            end

        else
            for i in 1:length(pivotElts) 
                pivotElt = pivotElts[i]
                diffPivotElt = pivotdiffElts[i]
                pivotColRow[pivotElt] = diffPivotElt
                dweight, dpath = diffPivotElt
                weight, path = pivotElt 
                entryTime[path] = max(weight, dweight)
            end
            #println("PCR: ",pivotColRow)
            for i in 1:length(basisElts)
                basisElt = basisElts[i]
                weight, path = basisElt 
                entryTime[path] = max(weight, 0)
            end
        end
        return entryTime, pivotColRow
    end

    function calculateEntryTime2(pivotElts, pivotdiffElts, basisElts, n, RP_1)
        #println("PIVOT ELTS: ", pivotElts)
        #println(length(pivotElts))
        entryTime = Dict()
        pivotColRow = Dict()
        if n == 1
            for i in 1:length(pivotElts) 
            
                pivotElt = pivotElts[i]
                diffPivotElt = pivotdiffElts[i]
                pivotColRow[pivotElt] = diffPivotElt
                dweight, dpath = diffPivotElt
                weight, path = pivotElt 
                entryTime[path] = max(weight, 0)
            end
            #println("PCR: ",pivotColRow)
            for i in 1:length(basisElts)
                basisElt = basisElts[i]
                weight, path = basisElt 
                entryTime[path] = max(weight, 0)
            end

        else
            for i in 1:length(pivotElts) 
                pivotElt = pivotElts[i]
                diffPivotElt = pivotdiffElts[i]
                pivotColRow[pivotElt] = diffPivotElt
                dweight, dpath = diffPivotElt
                weight, path = pivotElt 
                entryTime[path] = max(weight, dweight)
            end
            #println("PCR: ",pivotColRow)
            for i in 1:length(basisElts)
                basisElt = basisElts[i]
                boundaryRP = homology.dV2(basisElt[2],n+1,1)
                dWeight = 0 
                dPaths = collect(keys(boundaryRP))
                for p in dPaths 
                    pw = RP_1[p]
                    dWeight = max(dWeight, pw)
                end

                weight, path = basisElt 
                entryTime[path] = max(weight, dWeight)
            end
        end
        return entryTime, pivotColRow
    end

    function reverseIndicies(A, len)
        l = len + 1
        newA = []
        for a in A 
            push!(newA, l - a) 
        end
        return newA
    end

    function flipList(A)
        
    end
    #=
    ISSUE: NOT SELECTING BASIS ELEMENT ROWS CORRECTLY!!!
    =#
    function buildOmegaBases3(regularPaths,n)
        transformationMatrices = []
        reducedMatrices = []
        basisElts = [[] for i in 0:n]
        entryTimes = [Dict() for i in 0:n]
        dPivotElts = [[] for i in 0:n]
        pivotRowColRelation = [Dict() for i in 0:n]
        RP = regularPaths[2]
        RPList = sort(dictToList(RP))
        RP_1 = regularPaths[1]
        RP_1List = sort(dictToList(RP_1), rev=true)
        zeroPaths = []
        ETime = [Dict() for i in 0:n]
        ET = Dict()
        for vertex in RP_1List 
            push!(zeroPaths, [0, vertex])
            ET[vertex] = 0
        end

        entryTimes[1] = ET
        basisElts[1] = zeroPaths
        println(RP)
        M,wghtRng,eTime = buildMatrix(RP, RP_1, 1)
        ##println("eTime: ", eTime)
        ETime[1] = eTime
        Mreduced, T, nzcols, pivotcols, column_tracker,pivotR = columnReduceMatrix3(M)
        MreducedNEW, TN, nzcolsN, pivotcolsN, column_trackerN = columnReduceMatrix4(M)
        println("MreducedNEW: ", MreducedNEW)
        #println(column_tracker)
        #println("1-dim pivot: ", pivotcols)
        #println("1-dim: ", nzcols)
        pivotElts = getOmegaElts2(T,pivotcols,column_tracker,RPList)
        diffPivotElts = getPivotElts(Mreduced, pivotcols,RP_1List)
        dPivotElts[1] = diffPivotElts
        omegaElts = getOmegaElts2(T,nzcols,column_tracker,RPList)
        ET, pivotRC = calculateEntryTime(pivotElts, diffPivotElts, omegaElts,1)
        pivotRowColRelation[1]  = pivotRC
        entryTimes[2] = ET
        basisElts[2] = omegaElts
        println("basis ELTS: ", omegaElts)
        push!(transformationMatrices, T)
        push!(reducedMatrices, Mreduced)

        for i in 3:n
            #println("OMEGA ELTS: ", omegaElts)
            RP = regularPaths[i]
            RPList = sort(dictToList(RP))
            println(RPList)
            #println("RPLIST: ", RPList)
            RP_1 = regularPaths[i-1]
            RP_1List = sort(dictToList(RP_1), rev=true)
            if isempty(RPList) == false
                M,wghtRng,eTime = buildMatrix(RP, RP_1, i-1)
                ##println("eTime: ", eTime)
                ETime[i-1] = eTime

                ## SOMETHING IS WRONG WITH NZ COLS
                #println("ROW INFO")
                #println(column_tracker)
                #println(RP_1List)
                #println(nzcols)
                println("nzcols: ", nzcols)
                nzcols_tracked = column_tracker[nzcols]
                println("nzcols_tracked", nzcols_tracked)
                #println("nzcols_TRACKED: ",nzcols_tracked)
                rev_tracked = reverseIndicies(nzcols_tracked, length(column_tracker))
                println("rev_tracked: ", sort(rev_tracked))
                #println(rev_tracked)
                #println()
                println("M b4 reduce: ", M)
                rev_tracked = sort(rev_tracked)
                M = M[rev_tracked,:]
                Mc = deepcopy(M)
                println("M: ", M)
                Mreduced, T, nzcols, pivotcols,column_tracker,pivotR = columnReduceMatrix3(M)
                MreducedNEW, TN, nzcolsN, pivotcolsN, column_trackerN = columnReduceMatrix4(Mc)
                println("MreducedNEW: ", Mreduced)
                println("nz: ", findnz(Mreduced))
                println()
                #println("Column_TRACKER: ",column_tracker)
                #println("pivotcols: ", pivotcols)
                #println("nzcols: ", nzcols)
                pivotIndices = collect(keys(pivotR))
                PE =[]
                DIPE = [] 
                for index in pivotIndices 
                    higherElt = RPList[index]
                    push!(PE, higherElt)
                    lowerElt = RP_1List[pivotR[index]]
                    push!(DIPE, lowerElt)
                end

                pivotElts = getOmegaElts2(T,pivotcols,column_tracker,RPList)
                println("PIVOT ELTS: ",pivotElts)
                diffPivotElts = getPivotElts2(Mreduced, pivotcols,RP_1List, rev_tracked)
                println("DIFF PIVOT ELTS: ", diffPivotElts)
                #println("PIVOTCOLELTS: ", pivotElts)
                #println("DIFFPIVOTELTS: ", diffPivotElts)
                dPivotElts[i-1] = diffPivotElts
                omegaElts = getOmegaElts2(T,nzcols,column_tracker,RPList)
                #ET, pivotRC = calculateEntryTime(pivotElts, diffPivotElts, omegaElts,i)
                ET, pivotRC = calculateEntryTime(PE, DIPE, omegaElts,i)
                entryTimes[i] = ET
                #println("PIVOTRC: ", pivotRC)
                pivotRowColRelation[i-1] = pivotRC
                basisElts[i] = omegaElts
                push!(transformationMatrices, T)
                push!(reducedMatrices, Mreduced)
            end
        end
        #println("TOTOTOT: ", dPivotElts)
        return reducedMatrices, basisElts,entryTimes, dPivotElts, pivotRowColRelation,ETime
    end

    function persistencePathHomology(FG,n)
        filtration = shortestPathFiltration([FG[1],FG[2]], FG[3])

        allowedPaths = computeAllowTime([FG[1], FG[2]], FG[3], n)

        omegaMatrices = buildOmegaBases2(allowedPaths, n)

        persistenceHomology = [[] for i in 1:(length(FG[3])+1)]
        A1 =  sort(dictToList(allowedPaths[2]))
        wght = [p[1] for p in A1]
        indexRng = getRange(wght)
        #smithNormal = SNF.SNForm(omegaMatrices[1])
        numVertices = length(FG[1])
        push!(persistenceHomology[1], numVertices)
        for j in 1:length(indexRng)
            rng = indexRng[j]        
            smithNormalT = SNF.SNForm(omegaMatrices[1][:,1:rng[2]])

            diagonal_entries = []
            a,b = size(smithNormalT)
            for i in 1:min(a, b)
                x=  smithNormalT[i,i]
                if x == 0
                    break
                end
                push!(diagonal_entries, x)
            end
            H0 = numVertices - length(diagonal_entries)
            push!(persistenceHomology[j+1], H0)
        end
        SNFinfo = []
        for i in eachindex(omegaMatrices)
            omegaInfo = []
            allPaths =  sort(dictToList(allowedPaths[i+1]))
            wght = [p[1] for p in allPaths]
            indexRng = getRange(wght)
            for j in 1:length(indexRng)
                rng = indexRng[j]
                sf = homology.snfDiagonal(omegaMatrices[i][:,1:rng[2]])
                sf2 = SNF.SNForm(omegaMatrices[i])

                sf2 = sf2[:,1:rng[2]]
                a, b = size(sf2) 
                diagonal_entries = []
                for i in 1:min(a, b)
                    x= sf2[i,i]
                    if x == 0
                        break
                    end
                    push!(diagonal_entries, x)
                end


                torsion = homology.getTorsion(sf)
                push!(omegaInfo, [length(sf),length(diagonal_entries), torsion, size(omegaMatrices[i])[1]])                
            end
            push!(SNFinfo, omegaInfo)
        end
        return persistenceHomology, SNFinfo
    end

    function persistencePathHomology2(FG,n)
        filtration = shortestPathFiltration([FG[1],FG[2]], FG[3])

        allowedPaths = computeAllowTime([FG[1], FG[2]], FG[3], n)

        omegaMatrices = buildOmegaBases3(allowedPaths, n)

        persistenceHomology = [[] for i in 1:(length(FG[3])+1)]
        A1 =  sort(dictToList(allowedPaths[2]))
        wght = [p[1] for p in A1]
        indexRng = getRange(wght)
        #smithNormal = SNF.SNForm(omegaMatrices[1])
        numVertices = length(FG[1])
        push!(persistenceHomology[1], numVertices)
        for j in 1:length(indexRng)
            rng = indexRng[j]        
            smithNormalT = SNF.SNForm(omegaMatrices[1][:,1:rng[2]])

            diagonal_entries = []
            a,b = size(smithNormalT)
            for i in 1:min(a, b)
                x=  smithNormalT[i,i]
                if x == 0
                    break
                end
                push!(diagonal_entries, x)
            end
            H0 = numVertices - length(diagonal_entries)
            push!(persistenceHomology[j+1], H0)
        end
        SNFinfo = []
        for i in eachindex(omegaMatrices)
            omegaInfo = []
            allPaths =  sort(dictToList(allowedPaths[i+1]))
            wght = [p[1] for p in allPaths]
            indexRng = getRange(wght)
            for j in 1:length(indexRng)
                rng = indexRng[j]
                sf = homology.snfDiagonal(omegaMatrices[i][:,1:rng[1]])
                sf2 = SNF.SNForm(omegaMatrices[i])

                sf2 = sf2[:,1:rng[2]]
                a, b = size(sf2) 
                diagonal_entries = []
                for i in 1:min(a, b)
                    x= sf2[i,i]
                    if x == 0
                        break
                    end
                    push!(diagonal_entries, x)
                end


                torsion = homology.getTorsion(sf)
                push!(omegaInfo, [length(sf),length(diagonal_entries), torsion, size(omegaMatrices[i])[2], size(omegaMatrices[i])[1], rng[2]])                
            end
            push!(SNFinfo, omegaInfo)
        end
        return persistenceHomology, SNFinfo
    end

    function persistencePathHomology3(FG,n)
        filtration = shortestPathFiltration([FG[1],FG[2]], FG[3])

        allowedPaths = computeAllowTime([FG[1], FG[2]], FG[3], n)
        
        omegaMatrices,basisElts,entryTimes,dPivotElts, pivotRelation, ogETime = buildOmegaBases3(allowedPaths, n)
        Pers = [[] for i in 1:n]
        ###println("BasisElts: ", basisElts)
        ###println("EntryTimes: ",entryTimes)
        ###println(dPivotElts)
        ###println(pivotRelation)
        ##println("ogETime: ", ogETime)
        ##println(omegaMatrices)
        ##println("LENNNNN: ", length(omegaMatrices))
        for i in 1:(n-1)
            ogET = ogETime[i]
            ##println("HOMOLOGY DIM: ", i-1)
            listAP = collect(keys(allowedPaths[i+1]))
            listValues = collect(values(allowedPaths[i+1]))
            ET = entryTimes[i+1]
            pivotRel = pivotRelation[i]
            #println("ET: ", ET)
            #println("Basis Elts: ",basisElts[i+1])
            #println("dpivotElts: ",dPivotElts[i])
            #println("PIVOTREL: ", pivotRel)
            if isempty(basisElts[i+1]) == false
                for j in 1:length(listAP)
                    ap = listAP[j]
                    w = listValues[j]

                    if [w, ap] in basisElts[i+1]
                    else
                        if i == 1 
                            birth = 0
                            death = ET[ap]
                        else
                            #println("PATH: ", ap)
                            ##println("pivotRel: ", pivotRel)
                            et = max(pivotRel[[w, ap]][1],w)
                            #et = max(pivotRel[[w, ap]][1],w)
                            #birth = pivotRel[[w, ap]][1]

                            birth = entryTimes[i][pivotRel[[w, ap]][2]]
                            #birth = et
                            #death = ET[ap]
                            death = et
                            if i > 2
                                #println("og: ", ogETime[i-1][pivotRel[[w, ap]][2]])
                                birth = ogETime[i-1][pivotRel[[w, ap]][2]]
                            end
                            #println("birth, death, path: ", (birth, death,pivotRel[[w, ap]] ))
                        end

                        if birth < death 
                            #println("generates: ", ([w,ap]))
                            #println("FINITE BAR: birth and death; ", (birth,death))
                            push!(Pers[i], (birth, death))
                        end
                    end
                end
            end

            listAP_1 = collect(keys(allowedPaths[i]))
            listValues = collect(values(allowedPaths[i]))
            if isempty(basisElts[i]) == false
                for j in 1:length(listAP_1)
                    ap = listAP_1[j]
                    w = listValues[j] 
                    if i == 1 

                        if [0, [ap, ap]] in basisElts[i]
                            if [ap,ap] in dPivotElts[i]

                            else
                                push!(Pers[i], (0, Inf))
                            end
                        end
                    else
                        ET = entryTimes[i]
  
                        if [w, ap] in basisElts[i]
                            if [w, ap] in dPivotElts[i]

                            else
                                #println(dPivotElts[i])
                                #println(ET)
                                #println("[w,ap] generates: ", ([w,ap], ET[ap]))
                                # this line could be wrong!!
                                push!(Pers[i], (ET[ap], Inf))
                            end
                        end
                    end
                end
            end
            #println("PERS: ", Pers)
        end
        return Pers
    end

#=
    V = [1,2,3,4]
    E = [[1,2],[4,3],[2,3],[1,4]]
    W = [1,1,1,1]
    F = [1,2,4]
=#
#=
V = [1,2,3]
E = [[1,2],[2,3],[3,1]]
W = [1,1,1]
F = [1,2,3]
=#
#=
V = [1,2,3,4,5,6]
E = [[1,2],[3,2],[1,4],[3,4],[5,1],[5,2],[5,3],[5,4],[6,1],[6,2],[6,3],[6,4]]
W = [1,1,1,1,
    1,1,1,1,
    1,1,1,1]
F = [1,2,3]
=#

# Suspension of C3
    VSC3 = [1,2,3,4,5]
    ESC3 = [[1,2],[2,3],[3,1],[1,4],[2,4],[3,4],[1,5],[2,5],[3,5]]
    W = [1,1,1,1,1,1,1,1,1]
    F = [1,2,3]

    #=
    V = [1,2,3,4,5,6]
    E = [[1,2],[3,2],[3,4],[1,4],[6,1],[6,2],[6,3],[6,4],[5,1],[5,2],[5,3],[5,4]]
    W = [1 for e in E]
    F = [1,2,3]
    =#

    # Suspension of C3
    VSC3 = [1,2,3,4,5]
    ESC3 = [[1,2],[2,3],[3,1],[1,4],[2,4],[3,4],[1,5],[2,5],[3,5]]
    W = [1,1,1,1,1,1,1,1,1]
    F = [1,2,3]

    # Alternating 4-cycle
    VSC3 = [1,2,3,4]
    ESC3 = [[1,2],[3,2],[1,4],[3,4]]
    W = [1,1,1,1]
    F =[1,2,3]
  
    # C6 
    VSC3 = [1,2,3,4,5,6]
    ESC3 = [[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]]
    W = [1 for e in ESC3]
    F =[1,2,3,4,5,6]
 
    VSC3 = [1,2,3,4]
    ESC3 = [[1,2],[2,3],[3,4],[4,1]]
    W = [1 for e in ESC3]
    F =[1,2,3,4,5,6]

    VSC3 = [1,2,3,4,5]
    ESC3 = [[1,2],[2,3],[3,1],[1,4],[2,4],[3,4],[1,5],[2,5],[3,5]]
    W = [1,1,1,1,1,1,1,1,1]
    F = [1,2,3]

    VSC3 = [1,2,3,4]
    ESC3 = [[1,2],[2,4],[1,3],[3,4]]
    W = [1,2,3,4]
    F = [1,2,3,4,5,6]

    VSC3 = [1,2,3,4,5]
    ESC3 = [[1,2],[2,3],[3,1],[1,4],[2,4],[3,4],[1,5],[2,5],[3,5]]
    W = [1 for e in ESC3]
    F = [1,2,3,4,5,6]

    VSC3 = [1,2,3,4,5]
    ESC3 = [[1,2],[2,3],[3,4],[4,5],[5,1]]
    W = [1 for e in ESC3]
    F = [1,2,3,4,5,6]

    VSC3 = [1,2,3,4,5]
    ESC3 = [[1,2],[2,3],[3,4],[4,5],[5,1]]
    W = [1 for e in ESC3]
    F = [1,2,3,4,5,6]

    VSC3 = [1,2,3,4,5,6]
    ESC3 = [[1,2],[2,3],[3,1],[4,5],[5,6],[6,4]]
    W = [1 for e in ESC3]
    F = [1,2,3,4,5,6]

    wghtE = []
    for i in 1:length(ESC3)
        push!(wghtE, (W[i],ESC3[i]))
    end
    ##println(wghtE)

    filt = shortestPathFiltration([VSC3,wghtE],F)
    #println(filt)
    Ap = computeAllowTime([VSC3,wghtE],F,4)
    #println(Ap)
    #println(sort(dictToList(Ap[2])))
    #RP = regularPaths2([1,2,3,4],4)
    #updateRP = updateRegularPaths(Ap,RP, 4)
    #M,wghtRng = buildMatrix(Ap[2],Ap[1],1)

    #MP, T = (columnReduceMatrix(M, wghtRng))

   println(persistencePathHomology3([VSC3,wghtE, F], 4))


    ##println(sort(getOmegaBasis(T, updateRP[2])))
    ##println(MP)

    #=
        TM = buildOmegaBases2(Ap,4)

     for t in TM
        #println(t)
     end
for a in Ap 
        #println(a)
    end
    for f in filt
        #println(f)
    end

    sortedwghtE = sort(wghtE)
    #println(sortedwghtE)
    sortedE = []
    sortedW = []
    for i in 1:length(wghtE)
        push!(sortedE,sortedwghtE[i][2])
        push!(sortedW, sortedwghtE[i][1])
    end
    =#
end