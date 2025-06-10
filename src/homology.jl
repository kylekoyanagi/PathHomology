module homology
    using SparseArrays;
    #using LinearAlgebra;
    #using AbstractAlgebra;
    using BenchmarkTools;
    using Base.Threads;

    include("pathcomplex.jl")
    include("smith.jl")
    include("unionfind.jl")
    include("graph_preprocess.jl")
    using .SNF

    #==============================================================================================================
    getIndex input: 
    (1) A: an array 
    (2) a: an element to check if it in in a

    getIndex returns:
    (1) i: the position of the element a in a
    (2) -1: if the element a is not in A
    ==============================================================================================================#

    function getIndex(A,a)
        i = 0 
        for e in A
            i = i +1
            if (a == e) == true
                return i
            end
        end
        return -1
    end

    #==============================================================================================================
    getIndex input: 
    (1) A: an array 
    (2) a: an element to check if it in in a

    getIndex returns:
    (1) i: the position of the element a in a
    (2) -1: if the element a is not in A
    ==============================================================================================================#
    function getIndexV2(A,a)
        i = 0 
        for e in A
            i = i +1
            if issetequal(e,a) == true # checks if two sets are equal [i.e. [1,2] and [2,1] returns true]
                return i
            end
        end
        return -1
    end

    #==============================================================================================================
    # uniqueRows determines the unique paths in the differential and checks if each path is an allowed path
    uniqueRows input: 
    (1) A: allowed paths of length n-1
    (2) diff: The differential of allowed paths of length n 

    uniqueRows returns: 
    (1) uniqueRows: an array of the unique paths in the differential 
    (2) disallowed: an array of true/false valuse corresponding to true is the path is not allowed/false if the path is allowed                         
    ===============================================================================================================#
    function uniqueRows(A, diff)
        uniqueRows = []
        allowedPaths = []
        disallowed = []

        # turn paths in array form in to paths of string form [i.e. [1,0] t0 '10']
        for p in A
            allowedPaths = push!(allowedPaths, join(p))
        end

        for d in diff # iterate through each differential 
            for path in d # iterate through each path in the differential 

                if (join(path[2]) in uniqueRows) == false # check if the path has been added to uniqueRows

                    if (join(path[2]) in allowedPaths)  # check if the path is allowed
                        uniqueRows = push!(uniqueRows, join(path[2]))
                        disallowed = push!(disallowed, false)    
                    else
                        uniqueRows = push!(uniqueRows, join(path[2]))
                        disallowed = push!(disallowed, true)
                    end
                end
            end
        end
        return uniqueRows, disallowed
    end

    #==============================================================================================================
    # uniqueRows determines the unique paths in the differential and checks if each path is an allowed path
    uniqueRows input: 
    (1) A: allowed paths of length n-1
    (2) diff: The differential of allowed paths of length n 

    uniqueRows returns: 
    (1) uniqueRows: an array of the unique paths in the differential 
    (2) disallowed: an array of true/false valuse corresponding to true is the path is not allowed/false if the path is allowed                         
    ===============================================================================================================#
    function uniqueRowsV2(A, diff)
        uniqueRows = []
        allowedPaths = []
        disallowed = []

        # turn paths in array form in to paths of string form [i.e. [1,0] t0 '10']
        #for p in A
            #allowedPaths = push!(allowedPaths, join(p))
        #end
        pathKeys = keys(diff)
        for pathKey in pathKeys # iterate through each differential 
            diffKeys = keys(diff[pathKey])
            for diffKey in diffKeys # iterate through each path in the differential 
                if (diffKey in uniqueRows) == false # check if the path has been added to uniqueRows

                    if (diffKey in A)  # check if the path is allowed
                        uniqueRows = push!(uniqueRows, diffKey)
                        disallowed = push!(disallowed, false)    
                    else
                        uniqueRows = push!(uniqueRows, diffKey)
                        disallowed = push!(disallowed, true)
                    end
                end
            end
        end
        return uniqueRows, disallowed
    end

    #==============================================================================================================
    # isDegenerate determines if a path has a loop 
    isDegenerate input: 
    (1) path: a path in a differential 

    isDegenerate returns:
    (1) true: if the path is degenerate 
    (2) false: if the path is non-degenerate
    ===============================================================================================================#
    function isDegenerate(path)
        for i in 1:(length(path)-1)
            if path[i] == path[i+1]
                return true
            end
        end
        return false
    end

    #==============================================================================================================
    # Function used to print homology groups 
    ===============================================================================================================#
    function printSNF(snfH)
        println("The Homology Groups are: ")
        for i in 1:length(snfH)
            grp = ""
            freePart = snfH[i][1]
            torsionPart = snfH[i][2]
            if freePart > 0
                if freePart == 1
                    grp = grp*"Z"
                else
                    fp = string(freePart)
                    grp = grp*"Z^"*fp
                end
            end

            if isempty(torsionPart) == false
                tp = []
                for alpha in torsionPart
                    t = ""
                    a = string(alpha)
                    t = "Z/"*a
                    tp = push!(tp,t)
                end
                torsion = ""
                for j in 1:length(tp)
                    if j == length(tp)
                        torsion = torsion*tp[j]
                    else
                        torsion = torsion*tp[j]*"x"
                    end
                end
            end

            if (freePart == 0) && (isempty(torsionPart) == true)
                grp = grp*"0"
            elseif (freePart == 0) && (isempty(torsionPart) == false)
                grp = torsion
            elseif (freePart > 0) && (isempty(torsionPart) == false)
                grp = grp*"x"*torsion
            end
            print("H", i-1)
            println(" : ", grp)
        end
    end

    #==============================================================================================================
    # Calculates n-1 Allowed Paths
    A input:
    (1) X: a digraph 
    (2) n: maximum length of allowed paths to calculate n-1

    A returns: 
    (1) All allowed paths up to length n
    ==============================================================================================================#
    function A(X, n)
        A = []

        # Adds A_0 Paths to Allowed Paths A[1]
        sortedVertices = sort(string.(X.vertices))
        A = push!(A,sortedVertices)
        # Add A_1 Paths to Allowed Paths A[2]
        eList = []
        for e in X.edges
            eList = push!(eList, string.(e))
        end
        sortedEList = sort(eList)
        A = push!(A, sortedEList)
        if n > 1
            # Adds Allowable Paths A_2 to A_(n-1) (i.e. A[3] to A[n])
            for i in 3:n
                A_i = []
                for edge in A[2]
                    for path in A[i-1]
                        if edge[1] == path[length(path)] # if the start of an edge equals the end of the path
                            ap = vcat(path,[edge[2]])
                            A_i = push!(A_i,ap)
                        end
                    end
                end
                sortedA_i = sort(A_i)
                A = push!(A,sortedA_i) 
            end
        end
        return A
    end

    #==============================================================================================================
    # Calculate the differential of a path 
    d input: 
    (1) path: allowed path 
    (2) n: length of allowed path
    (3) pathSign: the sign of a path in a linear combination of paths

    d returns:
    (1) The differential for path 
    ==============================================================================================================#
    function dV2(path,n, pathSign)
        diff = Dict()
        if n == 1 # checks if path is a vertex ##check if vertex different condition!!!
            return diff[path] = 1
        else
            for i in 1:length(path) # calculate differential 
                sign = ((-1)^(i+1))*(pathSign)
                front = path[1:i-1]
                back = path[i+1:length(path)]
                together = vcat(front,back)
                isDegen = isDegenerate(together) # check if path is degenarate 
                if isDegen == true
                else
                    diff[together] = sign
                end
            end 
        end
            return diff
    end

    #==============================================================================================================
    # Calculates the differential for n
    d_n input: 
    (1) A: Allowed Paths array
    (2) n: nth differential to calculate

    d_n returns:
    (1) an array of for the nth differential
    ==============================================================================================================#
    function d_nV2(A,n)
        paths = A[n]
        differential = Dict()
        
        for path in paths
            diff = dV2(path,n,1)
            differential[path] = diff
        end
        return differential
    end

    function d_nV2P(A,n)
        paths = A[n]
        differential = Dict()
        if n == 1
            B = Channel{Tuple{String,Int64}}(length(paths))
            @threads for i in 1:length(paths)
                path = paths[i]
                diff = dV2(path,n,1)
                put!(B, (path, diff))
            end
            close(B)
            info = collect(B)  
        else
            B = Channel{Tuple{Vector{String},Dict{Any, Any}}}(length(paths))
            @threads for i in 1:length(paths)
                path = paths[i]
                diff = dV2(path,n,1)
                put!(B, (path, diff))
            end
            close(B)
            info = collect(B)
         end

        for (path,diff) in info 
            differential[path] = diff
        end

        return differential
    end
    #==============================================================================================================
    restrictedMatrixV2 input: 
    (1) A: An array of allowed paths [i.e. A(X,n)]
    (2) n: An index to specify which differentials to calculate [i.e. A[n]]

    restrictedMatrix returns:
    (1) a matrix representation of quotient C[n-1]/A[n-1]. An entry of +/-1 represents a path that is not in A[n-1].
    ===============================================================================================================#
    function restrictedMatrixV2(A,n)

        # calculate differential 
        diff = d_nV2P(A,n)

        # initialize matrix 
        rMatrix = []

        # get the unique paths of the differential and check if path is in A[n-1]
        uniRows, disallowed = uniqueRowsV2(A[n-1], diff)
        if n == 1
            rMatrix = []
            vertexList = keys(diff)
            for vertex in vertexList
                vertexCoeff = diff[vertex]
                rMatrix = push!(rMatrix, vertexCoeff)
            end
            return reduce(hcat, rMatrix), diff
        
        elseif isempty(diff) == true
            return rMatrix, diff
            
        else
            # iterate through each differential path 
            i = 1
            pathKeys = keys(diff)
            B = Channel{Vector{Int}}(length(pathKeys)) # output channel of images

            for pathKey in pathKeys
                col = zeros(length(uniRows))
                # iterate through each path of length n-1 in the differential 
                differential = diff[pathKey]
                diffKeys = keys(differential)
                for diffKey in diffKeys
                    substring = false
                    # iterate through each allowed path of length n-1
                    for apath in A[n-1]
                        for i in 1:length(apath)
                            if (apath[i] == diffKey[i])
                                substring = true
                            else
                                substring = false
                                break
                            end
                        end
                        if substring == true
                            break
                        end
                    end

                    # if path in the differential is not an allowed path add +/-1 depending on its sign 
                    if substring == true
                        
                    else 
                        index = getIndex(uniRows, diffKey) # get the position in the array that corresponds to the specific path
                        diffCoeff = differential[diffKey]
                        if diffCoeff == -1
                                col[index] = -1 
                        else
                                col[index] = 1
                        end
                    end
                end
                put!(B, col)
                #rMatrix = push!(rMatrix,col)
            end
        end

        close(B) # Close the channel when done

        rMatrix = collect(B) # Collect all images from the channel

        if length(rMatrix) == 1
            rMatrix = reshape(rMatrix[1],length(rMatrix[1]), 1)
            rMatrix = round.(Int64, rMatrix)
            return sparse(rMatrix), diff
        end
        rMatrix = reduce(hcat,rMatrix)
        rMatrix = round.(Int64, rMatrix)
        return sparse(rMatrix), diff
    end


    function restrictedMatrixV2P(A,n)

        # calculate differential 
        diff = d_nV2P(A,n)

        # initialize matrix 
        rMatrix = []

        # get the unique paths of the differential and check if path is in A[n-1]
        uniRows, disallowed = uniqueRowsV2(A[n-1], diff)
        if n == 1
            rMatrix = []
            vertexList = keys(diff)
            for vertex in vertexList
                vertexCoeff = diff[vertex]
                rMatrix = push!(rMatrix, vertexCoeff)
            end
            return reduce(hcat, rMatrix), diff
        
        elseif isempty(diff) == true
            return rMatrix, diff
            
        else
            # iterate through each differential path 
            i = 1
            pathKeys = collect(keys(diff))
            B = Channel{Vector{Int}}(length(pathKeys)) # output channel of images

            @threads for j in 1:length(pathKeys)
                pathKey = pathKeys[j]
                col = zeros(length(uniRows))
                # iterate through each path of length n-1 in the differential 
                differential = diff[pathKey]
                diffKeys = keys(differential)
                for diffKey in diffKeys
                    substring = false
                    # iterate through each allowed path of length n-1
                    for apath in A[n-1]
                        for i in 1:length(apath)
                            if (apath[i] == diffKey[i])
                                substring = true
                            else
                                substring = false
                                break
                            end
                        end
                        if substring == true
                            break
                        end
                    end

                    # if path in the differential is not an allowed path add +/-1 depending on its sign 
                    if substring == true
                        
                    else 
                        index = getIndex(uniRows, diffKey) # get the position in the array that corresponds to the specific path
                        diffCoeff = differential[diffKey]
                        if diffCoeff == -1
                                col[index] = -1 
                        else
                                col[index] = 1
                        end
                    end
                end
                put!(B, col)
                #rMatrix = push!(rMatrix,col)
            end
        end

        close(B) # Close the channel when done

        rMatrix = collect(B) # Collect all images from the channel

        if length(rMatrix) == 1
            rMatrix = reshape(rMatrix[1],length(rMatrix[1]), 1)
            rMatrix = round.(Int64, rMatrix)
            return sparse(rMatrix), diff
        end
        rMatrix = reduce(hcat,rMatrix)
        rMatrix = round.(Int64, rMatrix)
        return sparse(rMatrix), diff
    end
    #==============================================================================================================
    O_n input: 
    (1) A: An array of allowed paths [i.e. A(X,n)]
    (2) n: An index to specify which differentials to calculate [i.e. A[n]]

    restrictedMatrix returns:
    (1) O: An array of arrays where the elements in the array are dictionaries that represent elements in Omega n. 
    (2) OString: An array whose elements are dictionaries which give string representations of the paths in Omega n. (mainly used for indexing purposes and visualizing what is going on)
    ===============================================================================================================#
    function O_n(A,n)
        O = []  # an array to store the O_n's
        O_n = [] # an array to store all the elements in O_n  (i.e. stores comb)
        comb = Dict() # a dictionary to store an element in O_n
        OString = [] 
        oString = Dict()
        differentialList = [] # list to store the differentials 

        diff0 = d_nV2P(A,1)
        differentialList = push!(differentialList, diff0)

        diff1 = d_nV2P(A, 2)
        differentialList = push!(differentialList, diff1)

        # Add O_0 
        for v in A[1]    
            comb = Dict()
            comb[v] = 1
            push!(O_n, comb)

            # string representation 
            s = " + "*join(v)
            oString[v] =  s
        end 
        push!(OString, oString)
        push!(O, O_n)

        # Add O_1
        O_n = []
        oString = Dict()
        for e in A[2]
            comb = Dict()
            comb[e] = 1

            # string representation
            s = " + "*join(e)
            oString[e] = s
            push!(O_n, comb)
        end
        O = push!(O, O_n)
        push!(OString, oString)

        # Add O_n for O_2 to O_(n-1)
        for i in 3:n
            O_n = []
            oString = Dict()

            rMatrix, diff = restrictedMatrixV2(A,i) # Matrix representing Cn-1/An-1, the differentials for paths of length i
            push!(differentialList, diff) 
            pathKeys = collect(keys(diff)) # get a list of paths of length i, used for indexing later 
            if iszero(rMatrix) == false 
                rMatrix = sparse(Matrix{Int128}(rMatrix))
                rBasis = SNF.nullity(rMatrix)[2] # Calculate the matrix representing the nullspace of Cn-1/An-1

                # iterate through each column in the basis 
                for j in 1:size(rBasis)[2] 
                    comb = Dict() # dictionary to store elements in O_i 
                    v = rBasis[:,j] # column vector in the basis
                    linComb = [] # array to store the string representation of a path with its sign +/-
                    pathKey = [] # array to store the paths that the string represents 

                    # iterate through each row in the column
                    for k in 1:length(v) 
                        if v[k] != 0 
                            path = pathKeys[k]
                            comb[path] = v[k]
                            push!(pathKey, path)
                            sign = v[k]

                            # string representation 
                            if sign == 1
                                s = " + "*join(path)
                                push!(linComb, s)
                            else
                                s = " - "*join(path)
                                push!(linComb, s)
                            end
                        end
                    end 
                    lC = join(linComb)
                    oString[pathKey] = lC
                    push!(O_n, comb)               
                end
            else 
                for p in A[i]
                    comb = Dict()
                    comb[p] = 1
                    push!(O_n,comb)

                    s = " + "*join(p) 
                    oString[p] =  s
                end
            end
            push!(OString, oString)
            push!(O, O_n)
        end
        return O, OString, differentialList
    end

    function O_nP(A,n)
        O = []  # an array to store the O_n's
        O_n = [] # an array to store all the elements in O_n  (i.e. stores comb)
        comb = Dict() # a dictionary to store an element in O_n
        OString = [] 
        oString = Dict()
        differentialList = [] # list to store the differentials 

        diff0 = d_nV2P(A,1)
        differentialList = push!(differentialList, diff0)

        diff1 = d_nV2P(A, 2)
        differentialList = push!(differentialList, diff1)

        # Add O_0 
        for v in A[1]    
            comb = Dict()
            comb[v] = 1
            push!(O_n, comb)

            # string representation 
            s = " + "*join(v)
            oString[v] =  s
        end 
        push!(OString, oString)
        push!(O, O_n)

        # Add O_1
        O_n = []
        oString = Dict()
        for e in A[2]
            comb = Dict()
            comb[e] = 1

            # string representation
            s = " + "*join(e)
            oString[e] = s
            push!(O_n, comb)
        end
        O = push!(O, O_n)
        push!(OString, oString)

        # Add O_n for O_2 to O_(n-1)
        B = Channel{Tuple{Int64,Vector{Any},Dict{Any, Any},Dict{Any, Any}}}(n)
        @threads for i in 3:n
            O_n = []
            oString = Dict()
            rMatrix, diff = restrictedMatrixV2(A,i) # Matrix representing Cn-1/An-1, the differentials for paths of length i
            #push!(differentialList, diff) 
            pathKeys = collect(keys(diff)) # get a list of paths of length i, used for indexing later 
            if iszero(rMatrix) == false 
                rMatrix = sparse(Matrix{Int128}(rMatrix))
                rBasis = SNF.nullity(rMatrix)[2] # Calculate the matrix representing the nullspace of Cn-1/An-1

                # iterate through each column in the basis 
                for j in 1:size(rBasis)[2] 
                    comb = Dict() # dictionary to store elements in O_i 
                    v = rBasis[:,j] # column vector in the basis
                    linComb = [] # array to store the string representation of a path with its sign +/-
                    pathKey = [] # array to store the paths that the string represents 

                    # iterate through each row in the column
                    for k in 1:length(v) 
                        if v[k] != 0 
                            path = pathKeys[k]
                            comb[path] = v[k]
                            push!(pathKey, path)
                            sign = v[k]

                            # string representation 
                            if sign == 1
                                s = " + "*join(path)
                                push!(linComb, s)
                            else
                                s = " - "*join(path)
                                push!(linComb, s)
                            end
                        end
                    end 
                    lC = join(linComb)
                    oString[pathKey] = lC
                    push!(O_n, comb)               
                end
            else 
                for p in A[i]
                    comb = Dict()
                    comb[p] = 1
                    push!(O_n,comb)

                    s = " + "*join(p) 
                    oString[p] =  s
                end
            end
            #push!(OString, oString)
            #push!(O, O_n)


            put!(B, (i, O_n, oString, diff))
        end
        close(B)

        info = collect(B)
        sortedinfo = sort(info)
        for i in 1:length(sortedinfo)
            j, on, ostr, diff = sortedinfo[i]
            push!(OString, ostr)
            push!(O, on)
            push!(differentialList, diff)
        end

        return O, OString, differentialList
    end


    function pairDifferential(On,On_1, cancelTerms)
        OmegaList = collect(keys(On))
        Omega_1List = On_1
        #println(Omega_1List)
        omegaPairList = []
        omegaPairDict = Dict() 
        for omega in OmegaList
            newBoundary = Dict()
            omegaBoundary = On[omega]
            boundarySum = collect(keys(omegaBoundary))
            combination = []
            possibleCombs = []
            cancelTermsKeys = collect(keys(cancelTerms))

            if omega in cancelTermsKeys
                notAllowed = cancelTerms[omega]
            else 
                notAllowed = []
            end

            for elt in boundarySum
                added = false
                for omega_1 in Omega_1List
                    paths = collect(keys(omega_1))
                    signs = collect(values(omega_1))
                    if (elt in paths)

                        if (length(paths) == 1)
                            #index = findall(x->x==elt,paths)

                            if  signs[1]== omegaBoundary[elt]
                                newBoundary[elt] = 1
                                added = true
                                break
                            else
                                newBoundary[elt] = -1
                                added = true
                                break
                            end
                        else
                            push!(possibleCombs, (paths,signs))
                        end
                    end
                end

                if added == false
                    push!(combination, elt)
                end
            end


            loopPossibleCombs = copy(possibleCombs)
            for comb in loopPossibleCombs
                posPath = comb[1]
                signs = comb[2]
                len = 0
                pathComb = []
                for i in 1:length(posPath)
                    path = posPath[i]

                    if path in combination
                       len = len + 1 
                       push!(pathComb, (path, signs[i]))
                    end
                end

                if len == length(posPath)
                    sign = omegaBoundary[posPath[1]]
                    if sign == signs[1]
                        newBoundary[posPath] = 1
                        for p in posPath
                            filter!(x->x!=p,combination)
                        end
                        filter!(x->x!=comb, possibleCombs)

                    else
                        newBoundary[posPath] = -1
                        for p in posPath                           
                            filter!(x->x!=p,combination)
                        end
                        filter!(x->x!=comb, possibleCombs)
                    end
                end
            end
            loopPossibleCombs = copy(possibleCombs)
            loopCombination = copy(combination)
            for comb1 in loopCombination 
                for comb2 in loopCombination
                    if !(comb1 == comb2)
                        loopPossibleCombs = copy(possibleCombs)
                        comb1List = []
                        comb2List = []
                        for posComb in loopPossibleCombs
                            pths, sgns = posComb
                            if comb1 in pths
                                push!(comb1List, posComb)
                            end
                            if comb2 in pths
                                push!(comb2List, posComb)
                            end
                        end

                        for posComb1 in comb1List

                            pths1, sgns1 = posComb1
                            idx1 = getIndex(pths1, comb1)

                            removeElt1 = filter(x->x!=pths1[idx1], pths1)
                            #removeElt1sgns = filter(x->x!=sgns1[idx1], pths1)
                            #filter!(x->x!=(pths1[idx1], sgns1[idx1]), posComb1)
                            for posComb2 in comb2List

                                #idx2 = getIndex(pths2, comb2)
                                pths2, sgns2 = posComb2
                                #removeElt2 = filter(x->x!=(pths2[idx2], sgns2[idx2]),posComb2)
                                #println("removeElt, pths2", (removeElt1, pths2))
                                for remElt in removeElt1
                                    if remElt in pths2
                                    #if issubset(removeElt1,pths2) == true

                                        removeIdx = getIndex(pths2, removeElt1)

                                        sign1 = omegaBoundary[comb1]
                                        if sign1 == sgns1[idx1] 
                                            newBoundary[pths1] = 1
                                            for p in pths1
                                                filter!(x->x!=p,combination)
                                            end
                                            filter!(x->x!=comb1, possibleCombs)
                                        else
                                            newBoundary[pths1] = -1
                                            for p in pths1
                                                filter!(x->x!=p,combination)
                                            end
                                            filter!(x->x!=comb1, possibleCombs)
                                        end

                                        idx2 = getIndex(pths2, comb2)
                                        sign2 = omegaBoundary[comb2]
                                        if sign2 == sgns2[idx2]
                                            newBoundary[pths2] = 1
                                            for p in pths2
                                                filter!(x->x!=p,combination)
                                            end
                                            filter!(x->x!=comb2, possibleCombs)
                                        else
                                            newBoundary[pths2] = -1
                                            for p in pths2
                                                filter!(x->x!=p,combination)
                                            end
                                            filter!(x->x!=comb2, possibleCombs)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            #=
            loopPossibleCombs = copy(possibleCombs)
            for comb in loopPossibleCombs
                posPath = comb[1]
                signs = comb[2]
                len = 0
                pathComb = []
                println("posPath: ",posPath)
                cancelTerms = []
                contains = false
                for cancelTerm in notAllowed
                    cancelPath = cancelTerm[1]
                    cancelSign = cancelTerm[2]
                    println("CP: ",(cancelPath, posPath))
                    if cancelPath in posPath 
                        println("cancelPath: ",cancelPath)
                        push!(cancelTerms, cancelPath)
                        contains = true
                    end
                end
                if contains == true
                    for i in 1:length(posPath)
                        path = posPath[i]
                        if path in cancelTerms
                            #continue
                        else
                            if path in combination
                            len = len + 1 
                            push!(pathComb, (path, signs[i]))
                            end
                        end
                    end

                    if isempty(pathComb) == false
                        p,s = pathComb[1]
                        omegaComb = []
                        for om in p 
                            push!(omegaComb, om )
                        end
                        if s == omegaBoundary[p]
                            newBoundary[posPath] = 1
                            for p in posPath                         
                                filter!(x->x!=p,combination)
                            end
                            filter!(x->x!=comb, possibleCombs)
                            for cancelPath in cancelTerms
                                filter!(x->x!=(cancelPath, 1), notAllowed)
                            end
                        else
                            newBoundary[posPath] = -1
                            for p in posPath
                                filter!(x->x!=p,combination)
                            end
                            filter!(x->x!=comb, possibleCombs)
                            for cancelPath in cancelTerms
                                filter!(x->x!=(cancelPath, -1), notAllowed)
                            end
                        end
                    end
                end
            end
            =#

            omegaPairDict[omega] = newBoundary
            #push!(omegaPairList, newBoundary)
        end
        return omegaPairDict
    end

    #==============================================================================================================
    O_diffV3 input: 
    (1) On: A tuple from the function O_n. 
    (2) n: An index to specify which differentials to calculate [i.e. A[n]]

    restrictedMatrix returns:
    (1) OmegaDifferential: An array whose entries are the differential matrices for O_0 to O_(n-1)
    ===============================================================================================================#
    function O_diffV3(On, n)
        OInfo = On 
        OmegaDifferential = [] 

        # Add O_0
        O_0 = ones(1, length(OInfo[2][1]))
        O_0 = round.(Int64, O_0)
        #O_0 = matrix(ZZ, O_0)
        OmegaDifferential = push!(OmegaDifferential, O_0)
        
        emptyIndex = 0 # emptyIndex counts the number of non empty Omega's so we only loop through Omega's that have elements
        emptyCheck = OInfo[1]
        for arr in emptyCheck 
            if isempty(arr) == true
                break
            else
                emptyIndex = emptyIndex + 1
            end
        end

        for i in 1:(emptyIndex-1)
            OmegaD = []
            O = OInfo[1][i]
            OString = OInfo[2][i]
            On_1 = OInfo[1][i+1]
            On_1Dict = Dict()
            diffInfo = OInfo[3][i+1]
            notAllowed = Dict()
            # iterate through each element in Omega i+1
            for o in On_1 
                linearComb = collect(keys(o))
                differential = Dict()
                cancelTerms = []
                # iterate through each path in an element of Omega i+1
                for path in linearComb
                    pathSign = o[path] # sign of a path in an element in Omega

                    pathDiff = diffInfo[path] # differntial of path 
                    pathDiffKeys = keys(pathDiff) # paths in the differential of path 

                    # iterate through each path in the differential
                    for pathDiffKey in pathDiffKeys 
                        
                        differentialKeys = keys(differential) 
                        pathDiffCoeff = pathDiff[pathDiffKey] # coefficient of path in differential 

                        if pathDiffKey in differentialKeys # check if path has already been added into the differential 
                            currCoeff = differential[pathDiffKey]
                            coeff = pathSign*pathDiffCoeff
                            total = currCoeff + coeff 
                            differential[pathDiffKey] = total

                        else
                            coeff = pathDiffCoeff*pathSign
                            differential[pathDiffKey] = coeff
                        end
                    end
                end

                # iterate through each path in the differential dictionary
                for (path,sign) in differential 

                    if sign == 0 # remove all paths in the differential that have a 0 as its coefficient. 
                        #= 
                        IDEA: check if the sign is 0 in the above for loop. If it is remove it from the dict. 
                        =#
                        push!(cancelTerms, (path, 1))
                        push!(cancelTerms, (path, -1))
                        differential = delete!(differential, path)
                    end
                end
                notAllowed[linearComb] = cancelTerms
                On_1Dict[linearComb] = differential
                OmegaD = push!(OmegaD, differential)
            end
            OmegaMatrix = [] 
            OmegaDict = Dict()
            #println("On_1Dict: ", On_1Dict)

            if !(i == 1)
                pairedOmegaDict = pairDifferential(On_1Dict, O, notAllowed)
                #println("parieddict: ", pairedOmegaDict)
                #println(pairedOmegaDict)
                test = collect(keys(pairedOmegaDict))

                omega_elts = [k for k in test]

                for omega_elt in omega_elts
                    On_1String = collect(keys(OString))
                    OmegaNCol = zeros(length(keys(O)))
                    boundary_elts = collect(keys(pairedOmegaDict[omega_elt]))
                    
                    og_boundary = collect(keys(On_1Dict[omega_elt]))
                    goodBoundary_elts = []
                    matchings = Dict()
                    trackBoundary_elts = copy(boundary_elts)
                    for boundary_elt in boundary_elts
                        if typeof(boundary_elt) == Vector{String}
                            push!(goodBoundary_elts, boundary_elt)
                        else 
                            og_paths  = [] 
                            for path in boundary_elt 
                                if path in og_boundary
                                    push!(og_paths, path)
                                end
                            end
                            if length(boundary_elt) == length(og_paths)
                                push!(goodBoundary_elts, boundary_elt)

                            else
                                badPaths = []
                                for path in boundary_elt 
                                    if path in og_paths

                                    else
                                        push!(badPaths, path)
                                    end
                                end
                                matchings[boundary_elt] = (badPaths, length(badPaths))
                                for elt in trackBoundary_elts
                                    if !(elt == boundary_elt)
                                        if issubset(badPaths, elt) == true
                                            push!(goodBoundary_elts, boundary_elt)
                                            #push!(goodBoundary_elts, elt)
                                        end
                                    end
                                end
                            end
                        end
                    end

                    goodBoundary_elts = unique(goodBoundary_elts)

                    if isempty(matchings) == true

                    else
                        sortMatchings = sort(collect(matchings), by = x ->x.second[2], rev = true)
                        copySortMatchings = copy(sortMatchings)
                        for elt1 in copySortMatchings 
                            badPaths1, l1 = elt1[2]
                            allPaths1 = elt1[1]

                            hasMatch = []
                            contained = false
                            pathContained = Dict()
                            for path in badPaths1
                                pathContained[path] = []
                                for elt2 in copySortMatchings 
                                    badPaths2, l2 = elt2[2]
                                    allPaths2 = elt2[1]

                                    if !(allPaths1 == allPaths2)
                                        if path in badPaths2
                                            push!(pathContained[path], true)
                                        end
                                    end
                                end
                            end

                            for path in badPaths1 
                                if isempty(pathContained[path ]) == true
                                    filter!(x->x!=elt1, sortMatchings)   
                                end
                            end
                            #=
                            for elt2 in copySortMatchings 
                                badPaths2, l2 = elt2[2]
                                if !(badPaths1 == badPaths2)
                                    for path in badPaths1

                                        if path in badPaths2
                                            push!(hasMatch, 0)
                                        else
                                            push!(hasMatch, 1)
                                        end
                                    end
                                end
                            end
                            if length(hasMatch) == (length(copySortMatchings)-1)
                                println("elt1 remove: ", elt1)
                                filter!(x->x!=elt1, sortMatchings)   
                            end
                            =#
                        end
                        #=
                        if i == 4
                            sum = 0 
                            for sm in sortMatchings 
                                p,s = sm[2]
                                sum = sum + s
                            end
                            if sum > 2
                                println("omega: ", omega_elt)
                                println("og boundary: ", og_boundary)
                                println("sortMatchings: ", sortMatchings)
                                #println("mathcings: ", matchings)
                            end
                        end
                        =#
                        
                        good_omegas = [] 
                        for elt in sortMatchings 
                            push!(good_omegas, elt[1])
                        end

                        matchKeys = collect(keys(matchings))
                        matchingsMinusOmegas = setdiff(matchKeys, good_omegas)
                        for good_elt in good_omegas 
                            if good_elt in goodBoundary_elts 

                            else
                                push!(goodBoundary_elts, good_elt)
                            end
                        end

                        for bad_elt in matchingsMinusOmegas 
                            if bad_elt in goodBoundary_elts 
                                filter!(x->x!=bad_elt, goodBoundary_elts)
                            end
                        end
                        #temp = setdiff(goodBoundary_elts, matchingsMinusOmegas)
                        ## CURRENT ASSUMPTION: once you remove the omega elements that have a path that is not in any other omega element you are good to build the matrix
                        #goodBoundary_elts = temp 

                    end
                    goodBoundary_elts = unique(goodBoundary_elts)

                    #println("On_1Dict", On_1Dict)
                    for boundary_elt in goodBoundary_elts
                    #for boundary_elt in boundary_elts
                        if  typeof(boundary_elt) == Vector{String}
                            if length(boundary_elt) > 2
                                index = getIndex(On_1String, [boundary_elt])
                                if index == -1
                                    index = getIndex(On_1String, boundary_elt)
                                end

                            else 
                                index = getIndex(On_1String, boundary_elt)
                            end
                        else
                            index = getIndexV2(On_1String, boundary_elt) 
                        end
                        coeff = OmegaNCol[index]
                        OmegaNCol[index] = coeff + pairedOmegaDict[omega_elt][boundary_elt]
                    end
                    #=
                    if i == 4
                            println("omega_elt: ", omega_elt)
                            #println(On_1[omega_elt])
                            println("good boundary: ",goodBoundary_elts)
                            #println(boundary_elts)
                            println(OmegaNCol)
                            s = 0
                            for e in OmegaNCol
                                s = s + abs(e)
                            end
                            println((length(goodBoundary_elts), s))
                            println()
                    end
                    =#
                    unique_Dict = Dict()
                    for boundary_elt1 in goodBoundary_elts 
                        contained = false
                        unique_elts = copy(boundary_elt1)
                        for boundary_elt2 in goodBoundary_elts 
                            if !(boundary_elt1 == boundary_elt2)
                                # if first type is [1,2,3]
                                if typeof(boundary_elt1) == Vector{String} 
                                
                                    #=

                                    # if both are type [1,2,3], [1,3,4]
                                    if typeof(boundary_elt2) == Vector{String}
                                        if boundary_elt1 == boundary_elt2 
                                            filter!(x->x!=boundary_elt2 ,unique_elts)
                                            contained = true
                                            break

                                        else  

                                        end
                                    # type [1,2,3] and type [[1,2,3],[1,2,4], [1,4,5]]
                                    else
                                        if boundary_elt1 in boundary_elt2 
                                            filter
                                            contained = true 
                                            break
                                        end
                                    end
                                =#
                                # first type is [[1,2,3], [[1,3,4], [1,4,5]]]
                                else
                                    #second type is [1,2,3]
                                    if typeof(boundary_elt2) == Vector{String}
                                        if boundary_elt2 in boundary_elt1 
                                            filter!(x->x!=boundary_elt2, unique_elts)
                                            contained = true 
                                        end
                                    else
                                        for elt in boundary_elt1 
                                            if elt in boundary_elt2
                                                filter!(x->x!=elt, unique_elts)
                                                contained = true 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        unique_Dict[boundary_elt1] = unique_elts
                    end


                    omegaKeys = collect(keys(unique_Dict))
                    og_boundary = On_1Dict[omega_elt]
                    #println(og_boundary)
                    OmegaNCol = zeros(length(keys(O)))
                    for elt in omegaKeys
                        omegaSign = 0 
                        if typeof(unique_Dict[elt]) == Vector{String} 
                            uniqueElt = unique_Dict[elt]
                        else 
                            uniqueElt = unique_Dict[elt][1]
                        end

                        boundarySign = og_boundary[uniqueElt]
                        for j in 1:length(O)
                            omegaElt = O[j]
                            omegaEltKeys = collect(keys(omegaElt))
                            if typeof(elt) == Vector{String}  
                                if issetequal(omegaEltKeys,[elt]) == true 

                                    omegaSign = omegaElt[uniqueElt]
                                    break
                                end

                            else
                                if issetequal(omegaEltKeys,elt) == true 
                                    omegaSign = omegaElt[uniqueElt]
                                    break
                                end
                            end

                        end

                        if  typeof(elt) == Vector{String}
                            if length(elt) > 2
                                index = getIndex(On_1String, [elt])
                                if index == -1
                                    index = getIndex(On_1String, elt)
                                end

                            else 
                                index = getIndex(On_1String, elt)
                            end
                        else
                            index = getIndexV2(On_1String, elt) 
                        end

                        if omegaSign == boundarySign 
                            coeff = OmegaNCol[index]
                            OmegaNCol[index] = coeff + 1
                        else 
                            coeff = OmegaNCol[index]
                            OmegaNCol[index] = coeff -1
                        end
                    end
                    #=
                    if i == 2
                        println("omega_elt: ", omega_elt)
                        #println(On_1[omega_elt])
                        println("good boundary: ",goodBoundary_elts)
                        #println(boundary_elts)
                        println(OmegaNCol)
                        s = 0
                        for e in OmegaNCol
                            s = s + abs(e)
                        end
                        println((length(goodBoundary_elts), s))
                      
                        println()
                end
                =#
                    OmegaMatrix = push!(OmegaMatrix, OmegaNCol)
                end
            else


            # iterate through each differential in Omega n 
                for (K, ODiff) in On_1Dict
                    OmegaNCol = zeros(length(keys(O)))
                    OmegaKeys = keys(ODiff)
                    OStringKeys = keys(OString)

                    # interate through each elment in Omega n-1
                    for on in O 
                        OKeyList = keys(on)
                        combs = 1
                        linearCombination = [] 
                        #=
                        if i == 3
                            println(OKeyList)
                        end =#
                        # iterate each path in an element in Omega
                        for OKey in OKeyList
                            if typeof(OKey) == String
                                if [OKey] in OmegaKeys
                                    OKey = [OKey]
                                    if length(OKeyList) == combs
                                        index = getIndex(OStringKeys, OKey[1])
                                        
                                        ODiffKeySign = ODiff[OKey] # Sign of Path in differential
                                        OmegaKeySign = on[OKey[1]]

                                        if ODiffKeySign == OmegaKeySign
                                            coeff = OmegaNCol[index]
                                            OmegaNCol[index] = coeff + 1
                                        else
                                            coeff = OmegaNCol[index]
                                            OmegaNCol[index] = coeff - 1
                                        end
                                        ODiff = delete!(ODiff, OKey)
                                    end
                                end
                            else

                                if OKey in OmegaKeys # check if path in in the differential
                                    if (length(OKeyList) == combs) && (combs == 1) # check if the element in Omega is a linear combination or not
                                        index = getIndex(OStringKeys, OKey) 
                                        if index == -1 # checks if getIndex returns a valid index
                                            index = 1
                                            # if it doensn't find the valid index
                                            for OStringKey in OStringKeys
                                                if (typeof(OStringKey) == Vector{Any}) && (length(OStringKey) == 1)
                                                    Os = OStringKey[1]
                                                    if Os == OKey
                                                        break
                                                    else
                                                        index = index + 1
                                                    end
                                                else 
                                                    index = index + 1
                                                end
                                            end
                                        end

                                        ODiffKeySign = ODiff[OKey] # Sign of Path in differential
                                        OmegaKeySign = on[OKey] # Sign of path in Omega
                                        if ODiffKeySign == OmegaKeySign # if the two signs are equal 
                                            coeff = OmegaNCol[index] 
                                            OmegaNCol[index] = coeff + 1 # add 1
                                            #=
                                            if i == 3
                                                sum = 0 
                                                for p in OmegaNCol
                                                    sum = sum + abs(p)
                                                end
                                                println(OKeyList)
                                                println(sum)
                                            end=#
                                        else
                                            coeff = OmegaNCol[index]
                                            OmegaNCol[index] = coeff - 1 # subtract 1
                                           #= 
                                            if i == 4
                                                sum = 0 
                                                for p in OmegaNCol
                                                    sum = sum + abs(p)
                                                end
                                                println(OKeyList)
                                                println(sum)
                                            end
                                            =#
                                        end

                                        ODiff = delete!(ODiff, OKey) # remove key from differential 
                                        
                                    elseif length(OKeyList) == combs # If all the paths in the element in Omega are added
                                        linearCombination = push!(linearCombination, OKey)
                                        index = getIndexV2(OStringKeys, linearCombination) # get the index

                                        ODiffKeySign = ODiff[OKey] # Sign of Path in differential
                                        OmegaKeySign = on[OKey] # Sign of path in Omega

                                        if ODiffKeySign == OmegaKeySign # if the two signs are equal 
                                            coeff = OmegaNCol[index] 
                                            OmegaNCol[index] = coeff + 1 # add 1
                                            #=
                                            if i == 3
                                                sum = 0 
                                                for p in OmegaNCol
                                                    sum = sum + abs(p)
                                                end
                                                println(OKeyList)
                                                println(sum)
                                            end=#
                                        else
                                            coeff = OmegaNCol[index]
                                            OmegaNCol[index] = coeff - 1 # subtract 1
                                            #=
                                            if i == 3
                                                sum = 0 
                                                for p in OmegaNCol
                                                    sum = sum + abs(p)
                                                end
                                                println(OKeyList)
                                                println(sum)
                                            end=#
                                        end

                                        for p in linearCombination
                                            ODiff = delete!(ODiff, p)  # delete each path in the differential that was in the Omega Element
                                        end
                                    else
                                        linearCombination = push!(linearCombination, OKey)
                                        combs = combs + 1 
                                    end

                                # NEW PART TO SEE IF WORKS!!
                                else 
                                    linearCombination = push!(linearCombination, OKey)
                                    combs = combs + 1  
                                    if combs == (length(OKeyList) + 1)
                                        index = getIndexV2(OStringKeys, linearCombination) # get the index
                                        for okey in linearCombination
                                            if okey in OmegaKeys
                                                ODiffKeySign = ODiff[okey] # Sign of Path in differential
                                                OmegaKeySign = on[okey] # Sign of path in Omega
            
                                                if ODiffKeySign == OmegaKeySign # if the two signs are equal 
                                                    coeff = OmegaNCol[index] 
                                                    OmegaNCol[index] = coeff + 1 # add 1
                                                    #=
                                                    if i == 3
                                                        sum = 0 
                                                        for p in OmegaNCol
                                                            sum = sum + abs(p)
                                                        end
                                                        println(OKeyList)
                                                        println(sum)
                                                    end=#
                                                else
                                                    coeff = OmegaNCol[index]
                                                    OmegaNCol[index] = coeff - 1 # subtract 1
                                                    #=
                                                    if i == 3
                                                        sum = 0 
                                                        for p in OmegaNCol
                                                            sum = sum + abs(p)
                                                        end
                                                        println(OKeyList)
                                                        println(sum)
                                                    end =#
                                                end
                                            end
                                        end

                                        for p in linearCombination
                                            ODiff = delete!(ODiff, p)  # delete each path in the differential that was in the Omega Element
                                        end
                                    end
                                end
                        
                            end
                            OmegaKeys = keys(ODiff) # Update the keys in the differential 
                        end

                    end

                    if OmegaNCol in OmegaMatrix
                        for index in 1:length(OmegaMatrix)
                            if OmegaNCol == OmegaMatrix[index]
                                OmegaNCol = OmegaNCol - OmegaMatrix[index]
                                break
                            end
                        end
                    end
                    
                    if -1*OmegaNCol in OmegaMatrix
                        for index in 1:length(OmegaMatrix)
                            if -1*OmegaNCol == OmegaMatrix[index]
                                OmegaNCol = OmegaNCol + OmegaMatrix[index]
                                break
                            end
                        end
                    end
                    #=
                    if i == 3
                        sum = 0 
                        for p in OmegaNCol
                            sum = sum + abs(p)
                        end
                        println(sum)
                        println(OmegaNCol)
                        println()
                    end
                    =#
                    OmegaDict[K] = OmegaNCol
                    OmegaMatrix = push!(OmegaMatrix, OmegaNCol)
                end
            end


            if length(OmegaMatrix) == 1
                OmegaMatrix = reshape(OmegaMatrix[1], length(OmegaMatrix[1]),1)
            else
                OmegaMatrix = reduce(hcat, OmegaMatrix)
            end
            OmegaMatrixZZ = round.(Int128, OmegaMatrix)
            #OmegaMatrixZZ = matrix(ZZ,OmegaMatrixZZ) # convert to ring of integers
            OmegaMatrixZZ = SparseMatrixCSC{Int128}(OmegaMatrixZZ)
            OmegaDifferential = push!(OmegaDifferential, OmegaMatrixZZ) 
        end

        for j in emptyIndex:(n-1) # add empty array for all empty Omegas
            OmegaDifferential = push!(OmegaDifferential, [])
        end
        
        return OmegaDifferential
    end


    function O_diffV4(On, n)
        OInfo = On 
        OmegaDifferential = [] 

        # Add O_0
        O_0 = ones(1, length(OInfo[2][1]))
        O_0 = round.(Int64, O_0)
        #O_0 = matrix(ZZ, O_0)
        OmegaDifferential = push!(OmegaDifferential, O_0)
        
        emptyIndex = 0 # emptyIndex counts the number of non empty Omega's so we only loop through Omega's that have elements
        emptyCheck = OInfo[1]
        for arr in emptyCheck 
            if isempty(arr) == true
                break
            else
                emptyIndex = emptyIndex + 1
            end
        end

        for i in 1:(emptyIndex-1)
            OmegaD = []
            O = OInfo[1][i]
            OString = OInfo[2][i]
            On_1 = OInfo[1][i+1]
            On_1Dict = Dict()
            diffInfo = OInfo[3][i+1]
            notAllowed = Dict()
            # iterate through each element in Omega i+1
            for o in On_1 
                linearComb = collect(keys(o))
                differential = Dict()
                cancelTerms = []
                # iterate through each path in an element of Omega i+1
                for path in linearComb
                    pathSign = o[path] # sign of a path in an element in Omega

                    pathDiff = diffInfo[path] # differntial of path 
                    pathDiffKeys = keys(pathDiff) # paths in the differential of path 

                    # iterate through each path in the differential
                    for pathDiffKey in pathDiffKeys 
                        
                        differentialKeys = keys(differential) 
                        pathDiffCoeff = pathDiff[pathDiffKey] # coefficient of path in differential 

                        if pathDiffKey in differentialKeys # check if path has already been added into the differential 
                            currCoeff = differential[pathDiffKey]
                            coeff = pathSign*pathDiffCoeff
                            total = currCoeff + coeff 
                            differential[pathDiffKey] = total

                        else
                            coeff = pathDiffCoeff*pathSign
                            differential[pathDiffKey] = coeff
                        end
                    end
                end

                # iterate through each path in the differential dictionary
                for (path,sign) in differential 

                    if sign == 0 # remove all paths in the differential that have a 0 as its coefficient. 
                        #= 
                        IDEA: check if the sign is 0 in the above for loop. If it is remove it from the dict. 
                        =#
                        push!(cancelTerms, (path, 1))
                        push!(cancelTerms, (path, -1))
                        differential = delete!(differential, path)
                    end
                end
                notAllowed[linearComb] = cancelTerms
                On_1Dict[linearComb] = differential
                OmegaD = push!(OmegaD, differential)
            end
            OmegaMatrix = [] 
            OmegaDict = Dict()

            if !(i == 1)
                pairedOmegaDict = pairDifferential(On_1Dict, O, notAllowed)
                test = collect(keys(pairedOmegaDict))

                omega_elts = [k for k in test]

                for omega_elt in omega_elts
                    On_1String = collect(keys(OString))
                    OmegaNCol = zeros(length(keys(O)))
                    boundary_elts = collect(keys(pairedOmegaDict[omega_elt]))
                    
                    og_boundary = collect(keys(On_1Dict[omega_elt]))
                    goodBoundary_elts = []
                    matchings = Dict()
                    trackBoundary_elts = copy(boundary_elts)
                    for boundary_elt in boundary_elts
                        if typeof(boundary_elt) == Vector{String}
                            push!(goodBoundary_elts, boundary_elt)
                        else 
                            og_paths  = [] 
                            for path in boundary_elt 
                                if path in og_boundary
                                    push!(og_paths, path)
                                end
                            end
                            if length(boundary_elt) == length(og_paths)
                                push!(goodBoundary_elts, boundary_elt)

                            else
                                badPaths = []
                                for path in boundary_elt 
                                    if path in og_paths

                                    else
                                        push!(badPaths, path)
                                    end
                                end
                                matchings[boundary_elt] = (badPaths, length(badPaths))
                                for elt in trackBoundary_elts
                                    if !(elt == boundary_elt)
                                        if issubset(badPaths, elt) == true
                                            push!(goodBoundary_elts, boundary_elt)
                                            #push!(goodBoundary_elts, elt)
                                        end
                                    end
                                end
                            end
                        end
                    end

                    goodBoundary_elts = unique(goodBoundary_elts)

                    if isempty(matchings) == true

                    else
                        sortMatchings = sort(collect(matchings), by = x ->x.second[2], rev = true)
                        copySortMatchings = copy(sortMatchings)
                        for elt1 in copySortMatchings 
                            badPaths1, l1 = elt1[2]
                            allPaths1 = elt1[1]

                            hasMatch = []
                            contained = false
                            pathContained = Dict()
                            for path in badPaths1
                                pathContained[path] = []
                                for elt2 in copySortMatchings 
                                    badPaths2, l2 = elt2[2]
                                    allPaths2 = elt2[1]

                                    if !(allPaths1 == allPaths2)
                                        if path in badPaths2
                                            push!(pathContained[path], true)
                                        end
                                    end
                                end
                            end

                            for path in badPaths1 
                                if isempty(pathContained[path ]) == true
                                    filter!(x->x!=elt1, sortMatchings)   
                                end
                            end
                        end
                        
                        good_omegas = [] 
                        for elt in sortMatchings 
                            push!(good_omegas, elt[1])
                        end

                        matchKeys = collect(keys(matchings))
                        matchingsMinusOmegas = setdiff(matchKeys, good_omegas)
                        for good_elt in good_omegas 
                            if good_elt in goodBoundary_elts 

                            else
                                push!(goodBoundary_elts, good_elt)
                            end
                        end

                        for bad_elt in matchingsMinusOmegas 
                            if bad_elt in goodBoundary_elts 
                                filter!(x->x!=bad_elt, goodBoundary_elts)
                            end
                        end
                        #temp = setdiff(goodBoundary_elts, matchingsMinusOmegas)
                        ## CURRENT ASSUMPTION: once you remove the omega elements that have a path that is not in any other omega element you are good to build the matrix
                        #goodBoundary_elts = temp 

                    end
                    goodBoundary_elts = unique(goodBoundary_elts)
  
                    unique_Dict = Dict()
                    for boundary_elt1 in goodBoundary_elts 
                        contained = false
                        unique_elts = copy(boundary_elt1)
                        for boundary_elt2 in goodBoundary_elts 
                            if !(boundary_elt1 == boundary_elt2)
                                # if first type is [1,2,3]
                                if typeof(boundary_elt1) == Vector{String} 

                                # first type is [[1,2,3], [[1,3,4], [1,4,5]]]
                                else
                                    #second type is [1,2,3]
                                    if typeof(boundary_elt2) == Vector{String}
                                        if boundary_elt2 in boundary_elt1 
                                            filter!(x->x!=boundary_elt2, unique_elts)
                                            contained = true 
                                        end
                                    else
                                        for elt in boundary_elt1 
                                            if elt in boundary_elt2
                                                filter!(x->x!=elt, unique_elts)
                                                contained = true 
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        unique_Dict[boundary_elt1] = unique_elts
                    end


                    omegaKeys = collect(keys(unique_Dict))
                    og_boundary = On_1Dict[omega_elt]
                    #println(og_boundary)
                    OmegaNCol = zeros(length(keys(O)))
                    for elt in omegaKeys
                        omegaSign = 0 
                        if typeof(unique_Dict[elt]) == Vector{String} 
                            uniqueElt = unique_Dict[elt]
                        else 
                            uniqueElt = unique_Dict[elt][1]
                        end

                        boundarySign = og_boundary[uniqueElt]
                        for j in 1:length(O)
                            omegaElt = O[j]
                            omegaEltKeys = collect(keys(omegaElt))
                            if typeof(elt) == Vector{String}  
                                if issetequal(omegaEltKeys,[elt]) == true 

                                    omegaSign = omegaElt[uniqueElt]
                                    break
                                end

                            else
                                if issetequal(omegaEltKeys,elt) == true 
                                    omegaSign = omegaElt[uniqueElt]
                                    break
                                end
                            end

                        end

                        if  typeof(elt) == Vector{String}
                            if length(elt) > 2
                                index = getIndex(On_1String, [elt])
                                if index == -1
                                    index = getIndex(On_1String, elt)
                                end

                            else 
                                index = getIndex(On_1String, elt)
                            end
                        else
                            index = getIndexV2(On_1String, elt) 
                        end

                        if omegaSign == boundarySign 
                            coeff = OmegaNCol[index]
                            OmegaNCol[index] = coeff + 1
                        else 
                            coeff = OmegaNCol[index]
                            OmegaNCol[index] = coeff -1
                        end
                    end

                    OmegaMatrix = push!(OmegaMatrix, OmegaNCol)
                end
            else


            # iterate through each differential in Omega n 
                for (K, ODiff) in On_1Dict
                    OmegaNCol = zeros(length(keys(O)))
                    OmegaKeys = keys(ODiff)
                    OStringKeys = keys(OString)

                    # interate through each elment in Omega n-1
                    for on in O 
                        OKeyList = keys(on)
                        combs = 1
                        linearCombination = [] 
           
                        # iterate each path in an element in Omega
                        for OKey in OKeyList
                            if typeof(OKey) == String
                                if [OKey] in OmegaKeys
                                    OKey = [OKey]
                                    if length(OKeyList) == combs
                                        index = getIndex(OStringKeys, OKey[1])
                                        
                                        ODiffKeySign = ODiff[OKey] # Sign of Path in differential
                                        OmegaKeySign = on[OKey[1]]

                                        if ODiffKeySign == OmegaKeySign
                                            coeff = OmegaNCol[index]
                                            OmegaNCol[index] = coeff + 1
                                        else
                                            coeff = OmegaNCol[index]
                                            OmegaNCol[index] = coeff - 1
                                        end
                                        ODiff = delete!(ODiff, OKey)
                                    end
                                end
                            else

                                if OKey in OmegaKeys # check if path in in the differential
                                    if (length(OKeyList) == combs) && (combs == 1) # check if the element in Omega is a linear combination or not
                                        index = getIndex(OStringKeys, OKey) 
                                        if index == -1 # checks if getIndex returns a valid index
                                            index = 1
                                            # if it doensn't find the valid index
                                            for OStringKey in OStringKeys
                                                if (typeof(OStringKey) == Vector{Any}) && (length(OStringKey) == 1)
                                                    Os = OStringKey[1]
                                                    if Os == OKey
                                                        break
                                                    else
                                                        index = index + 1
                                                    end
                                                else 
                                                    index = index + 1
                                                end
                                            end
                                        end

                                        ODiffKeySign = ODiff[OKey] # Sign of Path in differential
                                        OmegaKeySign = on[OKey] # Sign of path in Omega
                                        if ODiffKeySign == OmegaKeySign # if the two signs are equal 
                                            coeff = OmegaNCol[index] 
                                            OmegaNCol[index] = coeff + 1 # add 1
        
                                        else
                                            coeff = OmegaNCol[index]
                                            OmegaNCol[index] = coeff - 1 # subtract 1
                                        end

                                        ODiff = delete!(ODiff, OKey) # remove key from differential 
                                        
                                    elseif length(OKeyList) == combs # If all the paths in the element in Omega are added
                                        linearCombination = push!(linearCombination, OKey)
                                        index = getIndexV2(OStringKeys, linearCombination) # get the index

                                        ODiffKeySign = ODiff[OKey] # Sign of Path in differential
                                        OmegaKeySign = on[OKey] # Sign of path in Omega

                                        if ODiffKeySign == OmegaKeySign # if the two signs are equal 
                                            coeff = OmegaNCol[index] 
                                            OmegaNCol[index] = coeff + 1 # add 1

                                        else
                                            coeff = OmegaNCol[index]
                                            OmegaNCol[index] = coeff - 1 # subtract 1
                                        end

                                        for p in linearCombination
                                            ODiff = delete!(ODiff, p)  # delete each path in the differential that was in the Omega Element
                                        end
                                    else
                                        linearCombination = push!(linearCombination, OKey)
                                        combs = combs + 1 
                                    end

                                # NEW PART TO SEE IF WORKS!!
                                else 
                                    linearCombination = push!(linearCombination, OKey)
                                    combs = combs + 1  
                                    if combs == (length(OKeyList) + 1)
                                        index = getIndexV2(OStringKeys, linearCombination) # get the index
                                        for okey in linearCombination
                                            if okey in OmegaKeys
                                                ODiffKeySign = ODiff[okey] # Sign of Path in differential
                                                OmegaKeySign = on[okey] # Sign of path in Omega
            
                                                if ODiffKeySign == OmegaKeySign # if the two signs are equal 
                                                    coeff = OmegaNCol[index] 
                                                    OmegaNCol[index] = coeff + 1 # add 1

                                                else
                                                    coeff = OmegaNCol[index]
                                                    OmegaNCol[index] = coeff - 1 # subtract 1
                                                end
                                            end
                                        end

                                        for p in linearCombination
                                            ODiff = delete!(ODiff, p)  # delete each path in the differential that was in the Omega Element
                                        end
                                    end
                                end
                        
                            end
                            OmegaKeys = keys(ODiff) # Update the keys in the differential 
                        end

                    end

                    if OmegaNCol in OmegaMatrix
                        for index in 1:length(OmegaMatrix)
                            if OmegaNCol == OmegaMatrix[index]
                                OmegaNCol = OmegaNCol - OmegaMatrix[index]
                                break
                            end
                        end
                    end
                    
                    if -1*OmegaNCol in OmegaMatrix
                        for index in 1:length(OmegaMatrix)
                            if -1*OmegaNCol == OmegaMatrix[index]
                                OmegaNCol = OmegaNCol + OmegaMatrix[index]
                                break
                            end
                        end
                    end
                    OmegaDict[K] = OmegaNCol
                    OmegaMatrix = push!(OmegaMatrix, OmegaNCol)
                end
            end


            if length(OmegaMatrix) == 1
                OmegaMatrix = reshape(OmegaMatrix[1], length(OmegaMatrix[1]),1)
            else
                OmegaMatrix = reduce(hcat, OmegaMatrix)
            end
            OmegaMatrixZZ = round.(Int128, OmegaMatrix)
            #OmegaMatrixZZ = matrix(ZZ,OmegaMatrixZZ) # convert to ring of integers
            OmegaMatrixZZ = SparseMatrixCSC{Int128}(OmegaMatrixZZ)
            OmegaDifferential = push!(OmegaDifferential, OmegaMatrixZZ) 
        end

        for j in emptyIndex:(n-1) # add empty array for all empty Omegas
            OmegaDifferential = push!(OmegaDifferential, [])
        end
        
        return OmegaDifferential
    end
    #==============================================================================================================
    snfDiagonal input: 
    (1) M: Matrix 

    snfDiagonal Returns: 
    (1) an array of the non-zero diagonal elements in the Smith Normal Form
    ==============================================================================================================#
    function snfDiagonal(M)
        if isempty(M) == false
            sf  = SNF.SNForm(M)
            #sf, T = SNF.smith(M)
            #sf, T = SNF.smithP(M)
            a, b = size(sf) 

            diagonal_entries = []
            for i in 1:min(a, b)
                x= sf[i,i]
                if x == 0
                    break
                end
                push!(diagonal_entries, x)
            end
            return diagonal_entries
        end
        return []
    end
    
    function getTorsion(arr)
        torsion = []
        for e in arr
            if !(abs(e) == 1)
                push!(torsion, abs(e))
            end
        end
        return torsion
    end

    #==============================================================================================================
    pathHomologyV2 input: 
    (1) X: A digraph X with specified vertices V and edges E
    (2) n: An integer that specifies the path homology to calculate up to.

    pathHomology returns: 
    (1) an array of the path homology for H0 to H(n-2) 
    ===============================================================================================================#
    function pathHomologyV2(X,n)

        # calculate the allowed paths in a digraph X
        allowedPaths = A(X,n)
        
        # initialize an array for the differentials and homology
        sfHomology = []

        # Calculate Omega 
        On = O_n(allowedPaths,n)

        # Calculate the Differentials 
        ODiff = O_diffV4(On, n)
        # calculate the path homology for H0
        zeroMtx = round.(Int128, zeros(1,1))
        zeroMtx = SparseMatrixCSC{Int128}(zeroMtx)
        sn1 = snfDiagonal(zeroMtx)
        sn2 = snfDiagonal(ODiff[2])
        rowLength = size(ODiff[2])[1]
    
        snfH0 = rowLength - length(sn2) - length(sn1)
        torsion = getTorsion(sn2) 
        push!(sfHomology,[snfH0,torsion])
        #= UNION FIND to get connected components, not ready for general graph input, 
        preprocessing work needs to be done to ensure labelling works. 
        torsion = 0
        snfH0 = unionfind.connectedComponents(X.vertices,X.edges)
        sfHomology = push!(sfHomology, [snfH0, torsion])
        =#

        # calculate the path homology for H1 to H(n-2)
        for i in 2:(n-1)
            if (isempty(ODiff[i]) == false) 
                if isempty(ODiff[i+1]) == true

                    sf = snfDiagonal(ODiff[i])

                    freePart = size(ODiff[i])[2] - length(sf)

                    torsionPart = []
                else
                    sf1 = snfDiagonal(ODiff[i])
                    sf2 = snfDiagonal(ODiff[i+1])

                    freePart = size(ODiff[i+1])[1] - length(sf2) - length(sf1)
                    torsionPart = getTorsion(sf2) 
                end
                sfHomology = push!(sfHomology, [freePart,torsionPart])
            else
                sfHomology = push!(sfHomology, [0,[]])
            end
        end

        return sfHomology
    end 

    #==============================================================================================================
    Version 3 is still in progress... 
    ===============================================================================================================#
    function pathHomologyV3(X,n)
        # preprocess graph 
        X.vertices, X.edges = graph_preprocess.cleanGraph(X.vertices,X.edges)

        # calculate the allowed paths in a digraph X
        allowedPaths = A(X,n)
        
        # initialize an array for the differentials and homology
        sfHomology = []

        # Calculate Omega 
        On = O_n(allowedPaths,n)
        # Calculate the Differentials 
        ODiff = O_diffV3(On, n)
        # calculate the path homology for H0
        #=
        zeroMtx = round.(Int128, zeros(1,1))
        zeroMtx = SparseMatrixCSC{Int128}(zeroMtx)
        sn1 = snfDiagonal(zeroMtx)
        sn2 = snfDiagonal(ODiff[2])
        rowLength = size(ODiff[2])[1]

        snfH0 = rowLength - length(sn2) - length(sn1)
        torsion = getTorsion(sn2) 
        =#
        #=UNION FIND to get connected components, not ready for general graph input, 
        preprocessing work needs to be done to ensure labelling works. =#
        torsion = []
        snfH0 = unionfind.connectedComponents(X.vertices,X.edges)
        sfHomology = push!(sfHomology, [snfH0, torsion])
    

        # calculate the path homology for H1 to H(n-2)
        for i in 2:(n-1)
            if (isempty(ODiff[i]) == false) 
                if isempty(ODiff[i+1]) == true

                    sf = snfDiagonal(ODiff[i])

                    freePart = size(ODiff[i])[2] - length(sf)

                    torsionPart = []
                else
                    sf1 = snfDiagonal(ODiff[i])
                    sf2 = snfDiagonal(ODiff[i+1])
                    freePart = size(ODiff[i+1])[1] - length(sf2) - length(sf1)
                    torsionPart = getTorsion(sf2) 
                end
                sfHomology = push!(sfHomology, [freePart,torsionPart])
            else
                sfHomology = push!(sfHomology, [0,[]])
            end
        end

        return sfHomology
    end 

    #==============================================================================================================
    HyperPathHomology input: 
    (1) X: A hypergraph X with specified vertices V and edges E
    (2) q: density q; that is the sequence of p n+1 vertices such that any q consequivive vertices of p lie in some edge.
    (3) n: An integer that specifies the path homology to calculate up to.

    HyperPathHomology returns: 
    (1) an array of the path homology for H0 to H(n-2) 
    ===============================================================================================================#
    function HyperPathHomology(X,q,n)

        # calculate the allowed paths in a digraph X
        allowedPaths = pathcomplex.buildpathcomplexV2(X,q,n)
        
        # initialize an array for the differentials and homology
        sfHomology = []

        # Calculate Omega 
        On = O_n(allowedPaths,n)

        # Calculate the Differentials 
        ODiff = O_diffV3(On, n)
        # calculate the path homology for H0
        zeroMtx = round.(Int128, zeros(1,1))
        zeroMtx = SparseMatrixCSC{Int128}(zeroMtx)
        sn1 = snfDiagonal(zeroMtx)
        sn2 = snfDiagonal(ODiff[2])
        rowLength = size(ODiff[2])[1]
        snfH0 = rowLength - length(sn2) - length(sn1)
        torsion = getTorsion(sn2) 
        sfHomology = push!(sfHomology, [snfH0, torsion])
        # calculate the path homology for H1 to H(n-2)
        for i in 2:(n-1)
            if (isempty(ODiff[i]) == false) 
                if isempty(ODiff[i+1]) == true

                    sf = snfDiagonal(ODiff[i])

                    freePart = size(ODiff[i])[2] - length(sf)

                    torsionPart = []
                else
                    sf1 = snfDiagonal(ODiff[i])
                    sf2 = snfDiagonal(ODiff[i+1])
                    freePart = size(ODiff[i+1])[1] - length(sf2) - length(sf1)
                    torsionPart = getTorsion(sf2) 
                end
                sfHomology = push!(sfHomology, [freePart,torsionPart])
            else
                sfHomology = push!(sfHomology, [0,[]])
            end
        end
        printSNF(sfHomology)
        return sfHomology
    end 

    function buildAdjacencyMatrix(X)
    end

    function addEdge(adj, u, v)
        push!(adj[u], v)
        push!(adj[v], u)
        return adj
    end

    function buildAdjacencyList(X)
        adj = [[] for _ in 1:length(X.vertices)]
        E = X.edges 
        for e in E 
            adj = addEdge(adj, e[1], e[2])
        end
        return adj
    end

    function tripleList(X)
        V = X.vertices
        adjList = buildAdjacencyList(X)
        triple_list = []
        for i in 1:(length(V)-1)
            u = V[i]
            uAdj = adjList[u]
            for j in 2:length(V)
                if i == j 
                    break
                end
                v = V[j]
                vAdj = adjList[v]
                W = intersect(uAdj,vAdj)

                if isempty(W) == false
                    push!(triple_list, (u,v, W))
                end
            end
        end
        return triple_list
    end

    function type1(triple_list, X)
        type1List = []
        E = X.edges
        for t in triple_list
            u = t[1]
            v = t[2]
            W = t[3]
            W1 = []
            for w in W
                if isempty(intersect([[u, w]], E)) == true 
                    break
                else 
                    if isempty(intersect([[w,v]], E)) == true
                        break
                    else 
                        push!(W1, w)
                    end
                end
            end
            if isempty(W1) == false
                push!(type1List, (u,v, W1))
            end
        end
        return type1List
    end

    function type2(triple_list, X)
        type2List = []
        E = X.edges
        for t in triple_list
            u = t[1]
            v = t[2]
            W = t[3]
            W2 = []
            for w in W
                if isempty(intersect([[v, w]], E)) == true 
                    break
                else 
                    if isempty(intersect([[w,u]], E)) == true
                        break
                    else 
                        push!(W2, w)
                    end
                end
            end
            if isempty(W2) == false
                push!(type1List, (u,v, W2))
            end
        end
        return type2List
    end

    function type3(triple_list, X)
        type3List = []
        E = X.edges

        for t in triple_list
            u = t[1]
            v = t[2]
            W = t[3]
            W3i = []
            W3j = []
            for i in 1:length(W)
                wi = W[i]
                for j in 1:length(W)
                    wj = W[j]
                    if i == j 
                        break
                    end
                    if isempty(intersect([[wi,u]], E)) == false
                        if isempty(intersect([[wi,v]], E)) == false
                            if isempty(intersect([[u,wj]], E)) == false
                                if isempty(intersect([[v,wj]], E)) == false
                                    push!(W3i, wi)
                                    push!(W3j, wj)
                                end
                            end
                        end
                    end

                end

            end
            if isempty(W3i) == false
                if isempty(W3j) == false
                    push!(type3List, (u,v,W3i, W3j))
                end
            end
        end
        return type3List
    end



end