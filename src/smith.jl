module SNF
using SparseArrays;
using LinearAlgebra;
#using AbstractAlgebra;
using BenchmarkTools;
using Base.Threads;

#=
    Functions for Smith Normal Form 
    Most of these functions are optimized for sparse matrix operations
=#
# Column Subtraction
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

# Scale a column by a factor
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

# Column Addition *if factor1 != +/-1, then scalecol! is used to permenatley scale col1*
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

# Scale a row by a factor
function scalerow!(x,row,scale)
    nzrows, nzcols, nzvals = findnz(x)
    nzindex = findall(x->x==row, nzrows)

    for i in nzindex
        col = nzcols[i]
        x[row,col] = scale*x[row,col]
    end
end

# Row Addition *if factor1 != +/-1, then scalerow! is used to permenatley scale row1*
function rowaddition!(M, row1, row2, factor1=1, factor2 = 1)
    nzrows, nzcols, nzvals = findnz(M)
    nzindex = findall(x->x==row2, nzrows)
    scalerow!(M, row1, factor1)
    for i in nzindex
        colindex = nzcols[i]
        M[row1,colindex] = M[row1,colindex] + factor2*M[row2,colindex]
    end
end

# Row Subtraction
function rowsubtraction!(M, row1, row2, factor=1)
    nzrows, nzcols, nzvals = findnz(M)
    nzindex = findall(x->x==row2, nzrows)

    for i in nzindex
        colindex = nzcols[i]
        M[row1,colindex] = M[row1,colindex] - factor*M[row2,colindex]
    end
end

# Swap Columns
function swapcol!(x,i,j)
    for k in axes(x,1) 
      idata = x[k,i]
      x[k,i] = x[k,j]
      x[k,j] = idata
    end
end

# Swap Rows
function swaprow!(x,i,j)
    for k in axes(x,2) 
        idata = x[i,k]
        x[i,k] = x[j,k]
        x[j,k] = idata
    end
end

#=
    nonzerocol(x,t)
    INPUT: 
        (1) x: sparse mxn matrix
        (2) t: pivot row/column
    OUTPUT: 
        (1) e: if a nonzero column entry exists its column number is returned; otherwise -1 is returned
=#
function nonzerocol(x,t)
    indices = findnz(dropzeros(x))[2]
    defferedCols = []
    uniIndices = unique(indices)

    # looks for a column with a +/-1 first, other columns are deffered
    for e in uniIndices
        if (e >= t) 
            if (1 in x[:,e]) || (-1 in x[:,e])
            return e
            else
                push!(defferedCols, e)
            end
        end
    end

    # if no column has a +/-1, we look at the deffered columns
    for e in defferedCols
        if (e >= t) 
            return e
        end
    end
    return -1
end

#=
    nonzerorow(x,t)
    INPUT: 
        (1) x: sparse mxn matrix
        (2) t: pivot row/column
    OUTPUT: 
        (1) e: if a nonzero row entry exists its row number is returned; otherwise -1 is returned
=#
function nonzerorow(x, j, t)
    nonzeroRow = findnz(x[:,j])[1]
    defferedRows = []

    # looks for a row with a +/-1 first, other rows are deffered
    for e in nonzeroRow
        if e >= t
            if abs(x[e,j]) == 1
                return e
            else
                push!(defferedRows, e)
            end
        end
    end

    # if no row with a +/-1 exist, look at the deffered rows
    for e in defferedRows
        if e >= t
            return e
        end
    end
    return -1
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

#=
    Smith Normal Form
    INPUT: 
        (1) M: sparse mxn matrix
        (2) with_transform: true if you want to calculate the column trasformation matrix; false otherwise
    OUTPUT: 
        (1) M: sparse mxn matrix in smith normal form
        (2) T: sparse nxn column transformation matrix 
=#
function smith(M, with_transform = true)

    m, n = size(M)
    min = minimum([m,n])

    if with_transform == true
        T = sparse(Matrix(1I, n, n))
    end

    for t in 1:min
        M = dropzeros(M)
        # find nonzero column, priority on columns with a +/- 1
        nonzeroCol = nonzerocol(M, t)

        if nonzeroCol == -1
            continue
        end

        # swap columns t and the nonzero column
        if !(t == nonzeroCol)
            swapcol!(M, t, nonzeroCol)
            if with_transform == true
                swapcol!(T, t, nonzeroCol)
            end
        end
        M = dropzeros(M)

        # find nonzero row, priority on row with a +/- 1
        nonzeroRow = nonzerorow(M, t, t)
        if nonzeroRow == -1
            continue
        end
        
        # swap row t and the nonzero row
        if !(t == nonzeroRow)
            swaprow!(M, t, nonzeroRow)
        end
        M = dropzeros(M)

        # clear everything under the pivot i.e. perform row operations
        for i in (t+1):m
            pivot = M[t,t]
            ai = M[i, t]

            # if the entry M[i,t] is nonzero, we need to clear it 
            if !(ai == 0)

                #= 
                If the pivot is 1, perform row subtraction: 
                    M[i,:] = M[i,:] - M[i,t]*M[t,:]
                to clear the entry M[i,t].
                =#
                if pivot == 1
                    rowsubtraction!(M,i,t,ai)

                #= 
                If the pivot divides M[i,t], calculate pivot/M[i,t], and perform row subtraction:
                    M[i,:] = M[i,:] - (pivot/M[i,t])*M[t,:]
                to clear the entry M[i,t]
                =#
                elseif (mod(ai,pivot) == 0)
                    factor = ai / pivot
                    rowsubtraction!(M,i,t,factor)
                
                #=
                Otherwise, use extended euclidean algorithm to get d,x,y such that:
                    d = M[t,t]*x + M[i,t]*y
                Then, use row addition:
                    M[i,:] = y*M[i,:] + x*M[t,:]
                to make the entry M[i,t] = d
                Finally, calculate pivot/d, scale row M[i,:] by pivot/d, and perform row subtraction:
                    M[i,:] = M[i,:] - M[t,:]
                to clear the entry M[i,t]
                =#
                else
                    d, x, y = gcdExt(pivot, ai)
                    rowaddition!(M, i, t, y, x)
                    factor = pivot / d
                    scalerow!(M,i, factor)
                    rowsubtraction!(M,i,t,1)
                end
            end
            M = dropzeros(M)
        end

        # clear everything to the right of the pivot i.e. perform column operations.
        for j in (t+1):n
                pivot = M[t,t]
                aj = M[t, j]

            # if the entry M[i,t] is nonzero, we need to clear it 
            if !(aj == 0)

                #= 
                If the pivot is 1, perform column subtraction: 
                    M[:,i] = M[:,i] - M[t,i]*M[:,t]
                to clear the entry M[t,i].
                =#
                if pivot == 1
                    colsubtraction!(M,j,t,aj)

                    if with_transform == true
                        colsubtraction!(T,j,t,aj)
                    end

                #= 
                If the pivot divides M[i,t], calculate pivot/M[i,t], and perform column subtraction:
                    M[:,i] = M[:,i] - (pivot/M[t,i])*M[:,t]
                to clear the entry M[t,i]
                =#                   
                elseif (mod(aj,pivot) == 0)
                    factor = aj / pivot
                        colsubtraction!(M,j,t,factor)

                    if with_transform == true
                        colsubtraction!(T,j,t,factor)
                    end

                #=
                Otherwise, use extended euclidean algorithm to get d,x,y such that:
                    d = M[t,t]*x + M[t,i]*y
                Then, use column addition:
                    M[:,i] = y*M[:,i] + x*M[:,t]
                to make the entry M[t,i] = d
                Finally, calculate pivot/d, scale row M[:,i] by pivot/d, and perform column subtraction:
                    M[:,i] = M[:,i] - M[:,t]
                to clear the entry M[t,i]
                =#                    
                else
                    d, x, y = gcdExt(pivot, aj)
                    coladdition!(M,j,t,y,x)

                    if with_transform == true
                        coladdition!(T,j,t,y,x)

                    end
                    factor = pivot / d
                    scalecol!(M, j, factor)
                    colsubtraction!(M,j,t,1)

                    if with_transform == true
                        scalecol!(T, j, factor)
                        colsubtraction!(T,j,t,1)
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
    Functions for Parallel Smith Normal Form
=#

# Partitions Rows/Columns of a Matrix Uniformly
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

# Parallel Column Subtraction
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

# Parallel Column Scaling
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

# Parallel Column Addition
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

# Parallel Row Scaling
function scalerow!(x,row,scale)
    nzrows, nzcols, nzvals = findnz(x)
    nzindex = findall(x->x==row, nzrows)

    for i in nzindex
        col = nzcols[i]
        x[row,col] = scale*x[row,col]
    end
end

# Parallel Row Addition
function rowadditionp!(M, row1, row2, factor1=1, factor2 = 1)
    nzrows, nzcols, nzvals = findnz(M)
    nzindex = findall(x->x==row2, nzrows)
    scalerow!(M, row1, factor1)
    for i in nzindex
        colindex = nzcols[i]
        M[row1,colindex] = M[row1,colindex] + factor2*M[row2,colindex]
    end
end

# Parallel Row Subtraction
function rowsubtractionp!(M, row1, row2, factor=1)
    nzrows, nzcols, nzvals = findnz(M)
    nzindex = findall(x->x==row2, nzrows)

    for i in nzindex
        colindex = nzcols[i]
        M[row1,colindex] = M[row1,colindex] - factor*M[row2,colindex]
    end
end

#=
    Parallel Smith Normal Form
    INPUT: 
        (1) M: sparse mxn matrix
        (2) with_transform: true if you want to calculate the column trasformation matrix; false otherwise
    OUTPUT: 
        (1) M: sparse mxn matrix in smith normal form
        (2) T: sparse nxn column transformation matrix 
    Parallelization Idea: 
        (1) Partition the rows/columns of a matrix uniformly. All partitions contain the pivot, t. 
        (2) Each thread looks at a specific part of the partition and clears those rows/columns in its parition.
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

        # find nonzero column, priority on columns with a +/- 1
        nonzeroCol = nonzerocol(M, t)

        if nonzeroCol == -1
            continue
        end

        # swap columns t and the nonzero column
        if !(t == nonzeroCol)
            swapcol!(M, t, nonzeroCol)
            if with_transform == true
                swapcol!(T, t, nonzeroCol)
            end
        end

        # find nonzero row, priority on row with a +/- 1
        nonzeroRow = nonzerorow(M, t, t)
        if nonzeroRow == -1
            continue
        end
        
        # swap row t and the nonzero row
        if !(t == nonzeroRow)
            swaprow!(M, t, nonzeroRow)
        end


        rowPartition = distMatrix2(k,t,m)
        parSize = length(rowPartition)
        C = Channel{Tuple{SparseMatrixCSC{Int128, Int128}, Vector{Int64}}}(length(rowPartition))
        
        # loop through the partition
        @threads for l in 1:parSize
            indicies = rowPartition[l]
            mm = length(indicies)
            
            # clear everything under the pivot in the partition i.e. perform row operations
            for i in indicies
                if !(i == t)
                    pivot = M[t,t]
                    ai = M[i, t]

                    # If the entry M[i,t] is nonzero, we need to clear it 
                    if !(ai == 0)
                        #= 
                        If the pivot is 1, perform column subtraction: 
                            M[:,i] = M[:,i] - M[t,i]*M[:,t]
                        to clear the entry M[t,i].
                        =#
                        if pivot == 1
                            rowsubtraction!(M,i,t,ai)
       
                        #= 
                        If the pivot divides M[i,t], calculate pivot/M[i,t], and perform row subtraction:
                            M[i,:] = M[i,:] - (pivot/M[i,t])*M[t,:]
                        to clear the entry M[i,t]
                        =#
                        elseif (mod(ai,pivot) == 0)
                            factor = ai / pivot
                            rowsubtraction!(M,i,t,factor)

                        #=
                        Otherwise, use extended euclidean algorithm to get d,x,y such that:
                            d = M[t,t]*x + M[i,t]*y
                        Then, use row addition:
                            M[i,:] = y*M[i,:] + x*M[t,:]
                        to make the entry M[i,t] = d
                        Finally, calculate pivot/d, scale row M[i,:] by pivot/d, and perform row subtraction:
                            M[i,:] = M[i,:] - M[t,:]
                        to clear the entry M[i,t]
                        =#
                        else
                            d, x, y = gcdExt(pivot, ai)
                            rowaddition!(M, i, t, y, x)
                            factor = pivot / d
                            scalerow!(M,i, factor)
                            rowsubtraction!(M,i,t,1)
                        end
                    end
                end
            end
        end
        close(C)

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

        B = Channel{Tuple{SparseMatrixCSC{Int128, Int128}, SparseMatrixCSC{Int128, Int128}, Vector{Any}}}(length(colPartition))
        # loop through the partition
        @threads for l in 1:parSize
            indicies = colPartition[l]
            nn = length(indicies)

            # clear everything to the right of the pivot in the partition i.e. perform column operations
            for index in 1:length(indicies)
                j = indicies[index]
                if !(j == t)
                    pivot = M[t,t]
                    aj = M[t, j]
                    pivotcol = cols[l]
                    pivotcolT = colsT[l]

                    # If the entry M[t,j] is nonzero, we need to clear it 
                    if !(aj == 0)
                        #= 
                        If the pivot is 1, perform column subtraction: 
                            M[:,i] = M[:,i] - M[t,i]*M[:,t]
                        to clear the entry M[t,i].
                        =#      
                        if pivot == 1
                            colsubtractionp!(M,j,pivotcol,aj)

                            if with_transform == true
                                colsubtractionp!(T,j,pivotcolT,aj)
                            end

                        #= 
                        If the pivot divides M[i,t], calculate pivot/M[i,t], and perform column subtraction:
                            M[:,i] = M[:,i] - (pivot/M[t,i])*M[:,t]
                        to clear the entry M[t,i]
                        =#     
                        elseif (mod(aj,pivot) == 0)
                            factor = aj / pivot
                                colsubtractionp!(M,j,pivotcol,factor)

                            if with_transform == true
                                colsubtractionp!(T,j,pivotcolT,factor)
                            end

                        #=
                        Otherwise, use extended euclidean algorithm to get d,x,y such that:
                            d = M[t,t]*x + M[t,i]*y
                        Then, use column addition:
                            M[:,i] = y*M[:,i] + x*M[:,t]
                        to make the entry M[t,i] = d
                        Finally, calculate pivot/d, scale row M[:,i] by pivot/d, and perform column subtraction:
                            M[:,i] = M[:,i] - M[:,t]
                        to clear the entry M[t,i]
                        =#  
                        else
                            d, x, y = gcdExt(pivot, aj)
                            coladditionp!(M,j,pivotcol,y,x)

                            if with_transform == true
                                coladditionp!(T,j,pivotcolT,y,x)

                            end
                            factor = pivot / d
                            scalecol!(M, j, factor)
                            colsubtractionp!(M,j,pivotcol,1)

                            if with_transform == true
                                scalecol!(T, j, factor)
                                colsubtractionp!(T,j,pivotcolT,1)
                            end
                        end
                    end
                end
            end
        end
        close(B)

    end
    if with_transform == true
        return M,T
    end
    return M
end

#= 
    Function using Smith Normal Form to get a Nullspace Basis
    INPUT: 
        (1) M: sparse mxn matrix 
    OUTPUT: 
        (1) nullity: the dimension of the nullspase
        (2) nullspaceBasis: a matrix whose columns are vectors of the nullspace basis
=# 
function nullity(M)
    S,T = smithP(M)

    m,n = size(S)
    l = minimum([m,n])
    r = 0

    # Calculate Rank i.e. # of nonzero diagonal entries
    for i in 1:l
        if !(S[i,i] == 0)
            r = r + 1
        end
    end

    # nullity = # of columns - rank
    nullity = n - r
    index = r+1
    nullspaceBasis = T[:, index:n]
    return nullity, nullspaceBasis
end

end