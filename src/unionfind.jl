module unionfind
using BenchmarkTools

    function find(Parent, i)
        if Parent[i] != i 
            Parent[i] = find(Parent ,Parent[i])
        end
        return Parent[i]
    end

    function union(Parent, rank , x, y)
        s1 = find(Parent, x)
        s2 = find(Parent, y)

        if s1 != s2
            if rank[s1] < rank[s2]
                Parent[s1] = s2
            elseif rank[s1] > rank[s2]
                Parent[s2] = s1
            else
                Parent[s2] = s1
                rank[s1] += 1
            end
        end
        return Parent, rank
    end

    function DSU(n)
        Parent = collect(1:n)
        rank = ones(Int,n)
        return Parent, rank
    end

    function kruskals_mst(V,E)
        Parent, rank = DSU(length(V))
        count = 0 
        cost = 0
        for e in E
            x = e[1]
            y = e[2]

            if find(Parent, x) != find(Parent, y)
                Parent, rank = union(Parent, rank, x,y)
                count += 1
                cost += 1
                if count  == length(V) - 1
                    break 
                end
            end
        end
        return Parent
    end

    function connectedComponents(V,E)
        return length(unique(kruskals_mst(V,E)))
    end
end