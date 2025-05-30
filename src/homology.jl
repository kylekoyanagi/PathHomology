module homology
    using SparseArrays;
    #using LinearAlgebra;
    #using AbstractAlgebra;
    using BenchmarkTools;
    using Base.Threads;

    include("pathcomplex.jl")
    include("smith.jl")
    include("unionfind.jl")
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

    #==============================================================================================================
    restrictedMatrixV2 input: 
    (1) A: An array of allowed paths [i.e. A(X,n)]
    (2) n: An index to specify which differentials to calculate [i.e. A[n]]

    restrictedMatrix returns:
    (1) a matrix representation of quotient C[n-1]/A[n-1]. An entry of +/-1 represents a path that is not in A[n-1].
    ===============================================================================================================#
    function restrictedMatrixV2(A,n)

        # calculate differential 
        diff = d_nV2(A,n)

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

        diff0 = d_nV2(A,1)
        differentialList = push!(differentialList, diff0)

        diff1 = d_nV2(A, 2)
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

            # iterate through each element in Omega i+1
            for o in On_1 
                linearComb = keys(o) 
                differential = Dict()

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
                        differential = delete!(differential, path)
                    end
                end

                On_1Dict[linearComb] = differential
                OmegaD = push!(OmegaD, differential)
            end

            OmegaMatrix = [] 
            OmegaDict = Dict()

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
        sf  = SNF.SNForm(M)
        #f, T = SNF.smith(M)
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
        ODiff = O_diffV3(On, n)
        # calculate the path homology for H0
        zeroMtx = round.(Int128, zeros(1,1))
        zeroMtx = SparseMatrixCSC{Int128}(zeroMtx)
        sn1 = snfDiagonal(zeroMtx)
        sn2 = snfDiagonal(ODiff[2])
        rowLength = size(ODiff[2])[1]
        snfH0 = rowLength - length(sn2) - length(sn1)
        torsion = getTorsion(sn2) 

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