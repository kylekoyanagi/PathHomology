using SparseArrays;
using Base.Threads
using Random

Random.seed!(1234) 

# Parameters
m, n = 1000, 1000    # matrix dimensions
nnz = 2000 # number of non-zero entries
val_range = 1:5       # range for random integers

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
    nonzerocol(x,t)
    INPUT: 
        (1) x: sparse mxn matrix
        (2) t: pivot row/column
    OUTPUT: 
        (1) e: if a nonzero column entry exists its column number is returned; otherwise -1 is returned
=#
function nonzerocol(x,t)
    dropzeros!(x)
    indices = findnz(x)[2]
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

function createRSM(m,n,nnz,val_range)
    # Generate random row and column indices for the non-zero values
    row_indices = rand(1:m, nnz)
    col_indices = rand(1:n, nnz)
    values = rand(val_range, nnz)

    # Create sparse matrix
    return sparse(row_indices, col_indices, values, m, n)
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

M = createRSM(m,n,nnz,val_range)

function colsubtraction!(x, col1, col2, factor = 1)
    rows = copy(rowvals(x))
    nz = nzrange(x, col2)
    for i in nz
        row = rows[i]
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
    rows = copy(rowvals(x))
    scalecol!(x, col1, factor1)
    for i in nzrange(x, col2)
        row = rows[i]
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

function colsubtractionC!(x, col1, col2, factor = 1)
    rows = copy(rowvals(x))
    nz = nzrange(x, col2)
    for i in nz
        row = rows[i]
        col1[row] = col1[row] - factor*x[row,col2]
    end
end

function scalecolC!(col,scale)
    indices, vals = findnz(col)
    for i in eachindex(indices)
        index = indices[i]
        value = vals[i]
        col[index] = value*scale
    end
    return col
end

function coladditionC!(x, col1, col2, factor1 = 1, factor2 = 1)
    rows = copy(rowvals(x))
    scalecolC!(col1, factor1)
    for i in nzrange(x, col2)
        row = rows[i]
        col1[row] = col1[row] + factor2*x[row,col2]
    end
end

function scalerowC!(row,scale)
    indicies, vals = findnz(row)
    for i in eachindex(indicies)
        index = indicies[i]
        value = vals[i]
        row[index] = value*scale
    end
    return row
end

function rowadditionC!(M, row1, row2, factor1=1, factor2 = 1)
    indicies, vals = findnz(M[row2,:])
    scalerowC!(row1, factor1)

    for i in eachindex(indicies)
        index = indicies[i]
        value = vals[i]
        row1[index] = row1[index] + factor2*M[row2,index]
    end
    return row1
end

function rowsubtractionC!(M, row1, row2, factor=1)
    indices, vals = findnz(M[row2,:])

    for i in eachindex(indices)
        index = indices[i]
        value = vals[i]
        row1[index] = row1[index] - factor*M[row2,index]
    end
    return row1
end

function clearRow(M, row, pivot)
    nzcols = copy(findnz(M[row,:])[1])

    B = Channel{Tuple{SparseVector{Int64, Int64}, Int64}}(max(size(M)[1],size(M)[2]))
    @threads for i in eachindex(nzcols)
    #for i in eachindex(nzcols)

        index = nzcols[i]
        aj = M[row, index]
        colj = M[:,index]
        pivotval = M[row,pivot]
        if !(pivot == index)
            if pivotval == 1
                colsubtractionC!(M,colj,pivot,aj)

            elseif (mod(aj,pivotval) == 0)
                factor = aj / pivotval

                colsubtractionC!(M,colj,pivot,factor)
            else
                d, x, y = gcdExt(pivotval, aj)
                coladditionC!(M,colj,pivot,y,x)

                factor = pivotval / d
                scalecolC!(colj, factor)
                colsubtractionC!(M,colj,pivot,1)

            end

        end
        put!(B, (colj, index))
    end
    close(B)
    sparseV = collect(B)

    for (col, index) in sparseV
        M[:,index] = col
    end
    return M
end

function clearCol(M, col, pivot)
    nzrows = copy(findnz(M[:,col])[1])

    B = Channel{Tuple{SparseVector{Int64, Int64}, Int64}}(max(size(M)[1],size(M)[2]))
    @threads for i in eachindex(nzrows)
        index = nzrows[i]
        ai = M[index,col]
        rowi = M[index,:]
        pivotval = M[pivot, col]

        if !(pivot == index)
            if pivotval == 1
                rowsubtractionC!(M,rowi,pivot,ai)

            elseif (mod(ai,pivotval) == 0)
                factor = ai / pivotval

                rowsubtractionC!(M,rowi,pivot,factor)
            else
                d, x, y = gcdExt(pivotval, ai)
                rowadditionC!(M,rowi,pivot,y,x)

                factor = pivotval / d
                scalerowC!(rowi, factor)
                rowsubtractionC!(M,rowi,pivot,1)

            end
        end
        put!(B, (rowi, index))
    end
    close(B)

    sparseV = collect(B)
    for (row, index) in sparseV
        M[index,:] = row
    end
    return M
end

function SNForm(M)
    m, n = size(M)
    min = minimum([m,n])
    k = nthreads()
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

        M = clearRow(M,t,t)

        M = clearCol(M,t,t)
    end
    return M
end

