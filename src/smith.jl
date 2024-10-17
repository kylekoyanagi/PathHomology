module SNF
using SparseArrays;
using LinearAlgebra;
#using AbstractAlgebra;
using BenchmarkTools;
using Base.Threads;

#=
# Functions for Smith Normal Form 
# Most of these functions are optimized matrix operations for sparse matrices
=#
function colsubtraction!(x, col1, col2, factor = 1)
    rows = rowvals(x)
    vals = nonzeros(x)
    m, n = size(x)
    nz = nzrange(x, col2)
    for i in nz
        row = rows[i]
        val = vals[i]
        x[row, col1] = x[row, col1] - factor*x[row,col2]
    end
end

function scalecol!(x,col,scale)
    rows = rowvals(x)
    vals = nonzeros(x)
    m, n = size(x)

    for i in nzrange(x,col)
        row = rows[i]
        val = vals[i]
        x[row,col] = scale*x[row,col]
    end
end

function coladdition!(x, col1, col2, factor1 = 1, factor2 = 1)
    rows = rowvals(x)
    vals = nonzeros(x)
    m, n = size(x)
    scalecol!(x, col1, factor1)
    for i in nzrange(x, col2)
        row = rows[i]
        val = vals[i]
        x[row, col1] = x[row, col1] + factor2*x[row,col2]
    end
end

function scalerow!(x,row,scale)
    nzrows, nzcols, nzvals = findnz(x)
    nzindex = findall(x->x==row, nzrows)

    for i in nzindex
        col = nzcols[i]
        x[row,col] = scale*x[row,col]
    end
end

function rowaddition!(M, row1, row2, factor1=1, factor2 = 1)
    nzrows, nzcols, nzvals = findnz(M)
    nzindex = findall(x->x==row2, nzrows)
    scalerow!(M, row1, factor1)
    for i in nzindex
        colindex = nzcols[i]
        M[row1,colindex] = M[row1,colindex] + factor2*M[row2,colindex]
    end
end

function rowsubtraction!(M, row1, row2, factor=1)
    nzrows, nzcols, nzvals = findnz(M)
    nzindex = findall(x->x==row2, nzrows)

    for i in nzindex
        colindex = nzcols[i]
        M[row1,colindex] = M[row1,colindex] - factor*M[row2,colindex]
    end
end

function swapcol!(x,i,j)
    for k in axes(x,1) 
      idata = x[k,i]
      x[k,i] = x[k,j]
      x[k,j] = idata
    end
end

function swaprow!(x,i,j)
    for k in axes(x,2) 
        idata = x[i,k]
        x[i,k] = x[j,k]
        x[j,k] = idata
    end
end

function nonzerocol(x,t)
    indices = findnz(dropzeros(x))[2]
    defferedCols = []
    uniIndices = unique(indices)
    for e in uniIndices
        if (e >= t) 
            if (1 in x[:,e]) || (-1 in x[:,e])
            return e
            else
                push!(defferedCols, e)
            end
        end
    end

    for e in defferedCols
        if (e >= t) 
            return e
        end
    end
    return -1
end

function nonzerorow(x, j, t)
    nonzeroCol = findnz(x[:,j])[1]
    defferedRows = []
    for e in nonzeroCol
        if e >= t
            if abs(x[e,j]) == 1
                return e
            else
                push!(defferedRows, e)
            end
        end
    end

    for e in defferedRows
        if e >= t
            return e
        end
    end
    return -1
end

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

#=
# Smith Normal Form
=#
function smith(M, with_transform = true)
    m, n = size(M)
    min = minimum([m,n])
    if with_transform == true
        T = sparse(Matrix(1I, n, n))
    end

    for t in 1:min
        M = dropzeros(M)
        # find pivot 
        nonzeroCol = nonzerocol(M, t)

        if nonzeroCol == -1
            continue
        end

        if !(t == nonzeroCol)
            swapcol!(M, t, nonzeroCol)
            if with_transform == true
                swapcol!(T, t, nonzeroCol)
            end
        end
        M = dropzeros(M)

        nonzeroRow = nonzerorow(M, t, t)
        if nonzeroRow == -1
            continue
        end
        
        if !(t == nonzeroRow)
            swaprow!(M, t, nonzeroRow)
        end
        M = dropzeros(M)
    for i in (t+1):m
            pivot = M[t,t]
            ai = M[i, t]
            if !(ai == 0)
                if pivot == 1
                    rowsubtraction!(M,i,t,ai)
                elseif (mod(ai,pivot) == 0)
                    factor = ai / pivot
                    rowsubtraction!(M,i,t,factor)
                else
                    d, x, y = gcdExt(pivot, ai)
                    #display([d,x,y])
                    rowaddition!(M, i, t, y, x)
                    factor = pivot / d
                    scalerow!(M,i, factor)
                    rowsubtraction!(M,i,t,1)
                    #rowsubtraction!(M,t,i,factor)
                end
            end
            M = dropzeros(M)
        end

        for j in (t+1):n
                pivot = M[t,t]
                aj = M[t, j]

            if !(aj == 0)
                if pivot == 1
                    colsubtraction!(M,j,t,aj)
                    #M[:, j] = M[:, j] .- aj .* M[:,t]
                    if with_transform == true
                        colsubtraction!(T,j,t,aj)
                        #T[:, j] = T[:, j] .- aj .* T[:,t]
                    end
                elseif (mod(aj,pivot) == 0)
                    factor = aj / pivot
                        colsubtraction!(M,j,t,factor)

                    #M[:, j] = M[:, j] .- factor .* M[:, t]

                    if with_transform == true
                        #T[:, j] = T[:, j] .- factor .* T[:, t]
                        colsubtraction!(T,j,t,factor)
                    end
                else
                    d, x, y = gcdExt(pivot, aj)
                    #M[:, j] = y .* M[:,j] .+ x .* M[:,t]
                    coladdition!(M,j,t,y,x)

                    if with_transform == true
                        #T[:, j] = y .* T[:,j] .+ x .* T[:,t]
                        coladdition!(T,j,t,y,x)

                    end
                    factor = pivot / d
                    #M[:, t] = M[:, t] .- factor .* M[:,j]
                    scalecol!(M, j, factor)
                    colsubtraction!(M,j,t,1)

                    #colsubtraction!(M,t,j,factor)
                    #swapcols!(M, t, j)
                    if with_transform == true
                        #T[:, t] = T[:, t] .- factor .* T[:,j]
                        scalecol!(T, j, factor)
                        colsubtraction!(T,j,t,1)
                        #colsubtraction!(T,t,j,factor)

                    end
                end
            end
            M = dropzeros(M)
        end
    end
    if with_transform == true
        return M,T
    end
    return M
end

#=
# Functions for Parallel Smith Normal Form
=#
function distMatrix2(k,t,m)
    partition = [] 
    d = m - t
    dd = m - d
    for z in 1:(k)
        subset = []
        push!(subset, t)
        partitionSize = floor((m - t - z)/k)
        for c in 0:(partitionSize)
                push!(subset, Int64(z + dd + k*c))
        end
        push!(partition, subset)
    end
    return partition
end

function partitionMatrixR2(M, partition)
    matrixPartition = []
    for subset in partition
        push!(matrixPartition, (M[subset,:],subset))
    end
    return matrixPartition
end

function combineMatrixR2!(M,Mlist)
    for (Mbar,indicies) in Mlist
        for i in 1:length(indicies)
            M[indicies[i],:] = Mbar[i,:]
        end
    end
end

function colsubtractionp!(x, col1, pivotcol, factor = 1)
    rows = rowvals(x)
    vals = nonzeros(x)
    m, n = size(x)
    nz = findnz(pivotcol)[1]
    #println((m,n,length(pivotcol)))
    for i in 1:length(nz)
        row = nz[i]
        # Ensure the row is valid for the matrix dimensions
        if row <= m && row >= 1 && row <= length(pivotcol)
            # Perform the subtraction on column col1 of matrix x
            x[row, col1] = x[row, col1] - factor * pivotcol[row]
        else
            println("Warning: Row $row is out of bounds!")
        end
        #row = rows[i]
        #val = vals[i]
        #x[row, col1] = x[row, col1] - factor*pivotcol[row]
    end
end

function scalecol!(x,col,scale)
    rows = rowvals(x)
    vals = nonzeros(x)
    m, n = size(x)

    for i in nzrange(x,col)
        row = rows[i]
        val = vals[i]
        x[row,col] = scale*x[row,col]
    end
end

function coladditionp!(x, col1, pivotcol, factor1 = 1, factor2 = 1)
    rows = rowvals(x)
    vals = nonzeros(x)
    m, n = size(x)
    scalecol!(x, col1, factor1)
    nz = findnz(pivotcol)[1]
    for i in 1:length(nz)
        row = nz[i]
        #val = vals[i]
        x[row, col1] = x[row, col1] + factor2*pivotcol[row]
    end
end

function scalerow!(x,row,scale)
    nzrows, nzcols, nzvals = findnz(x)
    nzindex = findall(x->x==row, nzrows)

    for i in nzindex
        col = nzcols[i]
        x[row,col] = scale*x[row,col]
    end
end

function rowadditionp!(M, row1, row2, factor1=1, factor2 = 1)
    nzrows, nzcols, nzvals = findnz(M)
    nzindex = findall(x->x==row2, nzrows)
    scalerow!(M, row1, factor1)
    for i in nzindex
        colindex = nzcols[i]
        M[row1,colindex] = M[row1,colindex] + factor2*M[row2,colindex]
    end
end

function rowsubtractionp!(M, row1, row2, factor=1)
    nzrows, nzcols, nzvals = findnz(M)
    nzindex = findall(x->x==row2, nzrows)

    for i in nzindex
        colindex = nzcols[i]
        M[row1,colindex] = M[row1,colindex] - factor*M[row2,colindex]
    end
end

#=
# Parallel Smith Normal Form
=#
function smithP(M, with_transform = true)
    m, n = size(M)
    min = minimum([m,n])
    k = nthreads()
    if with_transform == true
        T = sparse(Matrix(1I, n, n))
    end
    for t in 1:min
        M = dropzeros(M)
        # find pivot 
        nonzeroCol = nonzerocol(M, t)

        if nonzeroCol == -1
            continue
        end

        if !(t == nonzeroCol)
            swapcol!(M, t, nonzeroCol)
            if with_transform == true
                swapcol!(T, t, nonzeroCol)
            end
        end
        #M = dropzeros(M)

        nonzeroRow = nonzerorow(M, t, t)
        if nonzeroRow == -1
            continue
        end
        
        if !(t == nonzeroRow)
            swaprow!(M, t, nonzeroRow)
        end
        #M = dropzeros(M)
        rowPartition = distMatrix2(k,t,m)
        parSize = length(rowPartition)
        #rowMatrices = partitionMatrixR2(M,rowPartition)
       C = Channel{Tuple{SparseMatrixCSC{Int128, Int128}, Vector{Int64}}}(length(rowPartition))
        @threads for l in 1:parSize
            indicies = rowPartition[l]
            mm = length(indicies)
            for i in indicies
                if !(i == t)
                    #i = indicies[index]
                    pivot = M[t,t]
                    
                    ai = M[i, t]
                    if !(ai == 0)
                        if pivot == 1
                            rowsubtraction!(M,i,t,ai)
                        elseif (mod(ai,pivot) == 0)
                            factor = ai / pivot
                            rowsubtraction!(M,i,t,factor)
                        else
                            d, x, y = gcdExt(pivot, ai)
                            #display([d,x,y])
                            rowaddition!(M, i, t, y, x)
                            factor = pivot / d
                            scalerow!(M,i, factor)
                            rowsubtraction!(M,i,t,1)
                            #rowsubtraction!(M,t,i,factor)
                        end
                    end
                    #M = dropzeros(M)
                end
            end
            #put!(C,indicies)
        end
        close(C)
        #info = collect(C)
        #combineMatrixR2!(M,info)

        colPartition = distMatrix2(k,t,n)
        parSize = length(colPartition)

        pivotcol = M[:,t]
        pivotcolT = T[:,t]
        cols = []
        colsT = []
        for i in 1:length(colPartition)
            push!(cols,pivotcol)
            push!(colsT,pivotcolT)
        end
        #colMatrices = partitionMatrixC(M,T,colPartition)
        B = Channel{Tuple{SparseMatrixCSC{Int128, Int128}, SparseMatrixCSC{Int128, Int128}, Vector{Any}}}(length(colPartition))
        @threads for l in 1:parSize
            indicies = colPartition[l]
            nn = length(indicies)

            for index in 1:length(indicies)
                j = indicies[index]
                if !(j == t)
                    #j = indicies[index]
                    pivot = M[t,t]
                    aj = M[t, j]
                    pivotcol = cols[l]
                    pivotcolT = colsT[l]
                    if !(aj == 0)
                        if pivot == 1
                            colsubtractionp!(M,j,pivotcol,aj)
                            #M[:, j] = M[:, j] .- aj .* M[:,t]
                            if with_transform == true
                                colsubtractionp!(T,j,pivotcolT,aj)
                                #T[:, j] = T[:, j] .- aj .* T[:,t]
                            end
                        elseif (mod(aj,pivot) == 0)
                            factor = aj / pivot
                                colsubtractionp!(M,j,pivotcol,factor)

                            #M[:, j] = M[:, j] .- factor .* M[:, t]

                            if with_transform == true
                                #T[:, j] = T[:, j] .- factor .* T[:, t]
                                colsubtractionp!(T,j,pivotcolT,factor)
                            end
                        else
                            d, x, y = gcdExt(pivot, aj)
                            #M[:, j] = y .* M[:,j] .+ x .* M[:,t]
                            coladditionp!(M,j,pivotcol,y,x)

                            if with_transform == true
                                #T[:, j] = y .* T[:,j] .+ x .* T[:,t]
                                coladditionp!(T,j,pivotcolT,y,x)

                            end
                            factor = pivot / d
                            #M[:, t] = M[:, t] .- factor .* M[:,j]
                            scalecol!(M, j, factor)
                            colsubtractionp!(M,j,pivotcol,1)

                            #colsubtraction!(M,t,j,factor)
                            #swapcols!(M, t, j)
                            if with_transform == true
                                #T[:, t] = T[:, t] .- factor .* T[:,j]
                                scalecol!(T, j, factor)
                                colsubtractionp!(T,j,pivotcolT,1)
                                #colsubtraction!(T,t,j,factor)
                            end
                        end
                    end
                    #M = dropzeros(M)
                end
            end
            #put!(B, (Mbar,Tbar, indicies))
        end
        close(B)
        #info = collect(B)
        #combineMatrixC!(M,T,info)

    end
    if with_transform == true
        return M,T
    end
    return M
end

#= 
# Function using Smith Normal Form to get a Nullspace Basis
=# 
function nullity(M)
    #S, U, T = parSNF(M)
    #println(size(M))
    S,T = smithP(M)
    #S,T = smith(M)

    m,n = size(S)
    l = minimum([m,n])
    r = 0
    for i in 1:l
        if !(S[i,i] == 0)
            r = r + 1
        end
    end
    null = n - r
    index = r+1
    nullspaceB = T[:, index:n]
    return null, nullspaceB
end

end