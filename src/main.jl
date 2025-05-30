#==============================================================================================================#
# Load Packages 
#==============================================================================================================#
using SparseArrays;
#using LinearAlgebra;
#using AbstractAlgebra;
using BenchmarkTools;
using Base.Threads;

#==============================================================================================================#
# Load Files
#==============================================================================================================#
include("smith.jl")
include("homology.jl")
include("pathcomplex.jl")

using .SNF
using .homology

#==============================================================================================================#
# Define digraph Data Structure
#==============================================================================================================#
mutable struct digraph
    vertices:: Array
    edges:: Array
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
#==============================================================================================================#
# EXAMPLES
#==============================================================================================================#
# C3 Example
VC3 = [1,2,3]
EC3 = [[1,2],[2,3],[3,1]]
GC3 = digraph(VC3,EC3)

# C3,1 Example
VC31 = [1, 2, 3, 4]
EC31 = [[1, 2], [2, 4], [4, 3], [1, 3]]
GC31 = digraph(VC31, EC31)

# C2,2 Example 
VC22 = [1,2,3,4]
EC22 =[[1,2],[2,4],[1,3],[3,4]]
GC22 = digraph(VC22,EC22)

#homology.pathHomologyV2(GC22,3)

# Undirected Square Example
VUC4 = [1,2,3,4]
EUC4 = [[1,2], [2,1], [2,3], [3,2], [3,4], [4,3], [4,1], [1,4]]
GUC4 = digraph(VUC4, EUC4)

## Example RP2
VRP2 = [1,2,3,4,5,
    "a","b","c","d","e","f","g","h","i","j",
    "k","l","m","o","p","n","q","r","s","t",
    "z"]
ERP2 = [
    [1,2],[2,3],[3,4],[4,5],[5,1],
    ["a","b"],["b","c"],["c","d"],["d","e"],["e","f"],["f","g"],["g","h"],["h","i"],["i","j"],["j","a"],
    ["k","l"],["l","m"],["m","o"],["o","p"],["p","n"],["n","q"],["q","r"],["r","s"],["s","t"],["t","k"],
    [1,"a"],[2,"b"],[3,"c"],[4,"d"],[5,"e"],[1,"f"],[2,"g"],[3,"h"],[4,"i"],[5,"j"],
    ["a","k"],["b","l"],["c","m"],["d","o"],["e","p"],["f","n"],["g","q"],["h","r"],["i","s"],["j","t"],
    ["k","z"],["l","z"],["m","z"],["o","z"],["p","z"],["n","z"],["q","z"],["r","z"],["s","z"],["t","z"]
    ]

# Add other direction of edges
VRP2 = string.(VRP2)
undirectedERP2 = []
for e in ERP2
    dir = [string.(e[1]), string.(e[2])]
    push!(undirectedERP2, dir)
    otherDir = [string.(e[2]), string.(e[1])]
    push!(undirectedERP2, otherDir)
end
GRP2 = digraph(VRP2,undirectedERP2)

#@time homology.pathHomologyV2(GRP2,4)
#println(homology.pathHomologyV2(GRP2,4))
# Example SC4 
VSC4 = ["a", 1,2,3,4,"b"]
ESC4 = [["a", 1], ["a", 2], ["a", 3], ["a", 4],
    [1, 2], [2, 3], [3, 4], [1, 4],
    ["b", 1], ["b", 2], ["b", 3], ["b", 4]]
GSC4 = digraph(VSC4,ESC4)

# 3-D Torus Example
V3DT = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125]
E3DT = Any[(1, 1), (1, 2), (1, 5), (1, 6), (1, 21), (1, 26), (1, 101), (2, 1), (2, 2), (2, 3), (2, 7), (2, 22), (2, 27), (2, 102), (3, 2), (3, 3), (3, 4), (3, 8), (3, 23), (3, 28), (3, 103), (4, 3), (4, 4), (4, 5), (4, 9), (4, 24), (4, 29), (4, 104), (5, 1), (5, 4), (5, 5), (5, 10), (5, 25), (5, 30), (5, 105), (6, 1), (6, 6), (6, 7), (6, 10), (6, 11), (6, 31), (6, 106), (7, 2), (7, 6), (7, 7), (7, 8), (7, 12), (7, 32), (7, 107), (8, 3), (8, 7), (8, 8), (8, 9), (8, 13), (8, 33), (8, 108), (9, 4), (9, 8), (9, 9), (9, 10), (9, 14), (9, 34), (9, 109), (10, 5), (10, 6), (10, 9), (10, 10), (10, 15), (10, 35), (10, 110), (11, 6), (11, 11), (11, 12), (11, 15), (11, 16), (11, 36), (11, 111), (12, 7), (12, 11), (12, 12), (12, 13), (12, 17), (12, 37), (12, 112), (13, 8), (13, 12), (13, 13), (13, 14), (13, 18), (13, 38), (13, 113), (14, 9), (14, 13), (14, 14), (14, 15), (14, 19), (14, 39), (14, 114), (15, 10), (15, 11), (15, 14), (15, 15), (15, 20), (15, 40), (15, 115), (16, 11), (16, 16), (16, 17), (16, 20), (16, 21), (16, 41), (16, 116), (17, 12), (17, 16), (17, 17), (17, 18), (17, 22), (17, 42), (17, 117), (18, 13), (18, 17), (18, 18), (18, 19), (18, 23), (18, 43), (18, 118), (19, 14), (19, 18), (19, 19), (19, 20), (19, 24), (19, 44), (19, 119), (20, 15), (20, 16), (20, 19), (20, 20), (20, 25), (20, 45), (20, 120), (21, 1), (21, 16), (21, 21), (21, 22), (21, 25), (21, 46), (21, 121), (22, 2), (22, 17), (22, 21), (22, 22), (22, 23), (22, 47), (22, 122), (23, 3), (23, 18), (23, 22), (23, 23), (23, 24), (23, 48), (23, 123), (24, 4), (24, 19), (24, 23), (24, 24), (24, 25), (24, 49), (24, 124), (25, 5), (25, 20), (25, 21), (25, 24), (25, 25), (25, 50), (25, 125), (26, 1), (26, 26), (26, 27), (26, 30), (26, 31), (26, 46), (26, 51), (27, 2), (27, 26), (27, 27), (27, 28), (27, 32), (27, 47), (27, 52), (28, 3), (28, 27), (28, 28), (28, 29), (28, 33), (28, 48), (28, 53), (29, 4), (29, 28), (29, 29), (29, 30), (29, 34), (29, 49), (29, 54), (30, 5), (30, 26), (30, 29), (30, 30), (30, 35), (30, 50), (30, 55), (31, 6), (31, 26), (31, 31), (31, 32), (31, 35), (31, 36), (31, 56), (32, 7), (32, 27), (32, 31), (32, 32), (32, 33), (32, 37), (32, 57), (33, 8), (33, 28), (33, 32), (33, 33), (33, 34), (33, 38), (33, 58), (34, 9), (34, 29), (34, 33), (34, 34), (34, 35), (34, 39), (34, 59), (35, 10), (35, 30), (35, 31), (35, 34), (35, 35), (35, 40), (35, 60), (36, 11), (36, 31), (36, 36), (36, 37), (36, 40), (36, 41), (36, 61), (37, 12), (37, 32), (37, 36), (37, 37), (37, 38), (37, 42), (37, 62), (38, 13), (38, 33), (38, 37), (38, 38), (38, 39), (38, 43), (38, 63), (39, 14), (39, 34), (39, 38), (39, 39), (39, 40), (39, 44), (39, 64), (40, 15), (40, 35), (40, 36), (40, 39), (40, 40), (40, 45), (40, 65), (41, 16), (41, 36), (41, 41), (41, 42), (41, 45), (41, 46), (41, 66), (42, 17), (42, 37), (42, 41), (42, 42), (42, 43), (42, 47), (42, 67), (43, 18), (43, 38), (43, 42), (43, 43), (43, 44), (43, 48), (43, 68), (44, 19), (44, 39), (44, 43), (44, 44), (44, 45), (44, 49), (44, 69), (45, 20), (45, 40), (45, 41), (45, 44), (45, 45), (45, 50), (45, 70), (46, 21), (46, 26), (46, 41), (46, 46), (46, 47), (46, 50), (46, 71), (47, 22), (47, 27), (47, 42), (47, 46), (47, 47), (47, 48), (47, 72), (48, 23), (48, 28), (48, 43), (48, 47), (48, 48), (48, 49), (48, 73), (49, 24), (49, 29), (49, 44), (49, 48), (49, 49), (49, 50), (49, 74), (50, 25), (50, 30), (50, 45), (50, 46), (50, 49), (50, 50), (50, 75), (51, 26), (51, 51), (51, 52), (51, 55), (51, 56), (51, 71), (51, 76), (52, 27), (52, 51), (52, 52), (52, 53), (52, 57), (52, 72), (52, 77), (53, 28), (53, 52), (53, 53), (53, 54), (53, 58), (53, 73), (53, 78), (54, 29), (54, 53), (54, 54), (54, 55), (54, 59), (54, 74), (54, 79), (55, 30), (55, 51), (55, 54), (55, 55), (55, 60), (55, 75), (55, 80), (56, 31), (56, 51), (56, 56), (56, 57), (56, 60), (56, 61), (56, 81), (57, 32), (57, 52), (57, 56), (57, 57), (57, 58), (57, 62), (57, 82), (58, 33), (58, 53), (58, 57), (58, 58), (58, 59), (58, 63), (58, 83), (59, 34), (59, 54), (59, 58), (59, 59), (59, 60), (59, 64), (59, 84), (60, 35), (60, 55), (60, 56), (60, 59), (60, 60), (60, 65), (60, 85), (61, 36), (61, 56), (61, 61), (61, 62), (61, 65), (61, 66), (61, 86), (62, 37), (62, 57), (62, 61), (62, 62), (62, 63), (62, 67), (62, 87), (63, 38), (63, 58), (63, 62), (63, 63), (63, 64), (63, 68), (63, 88), (64, 39), (64, 59), (64, 63), (64, 64), (64, 65), (64, 69), (64, 89), (65, 40), (65, 60), (65, 61), (65, 64), (65, 65), (65, 70), (65, 90), (66, 41), (66, 61), (66, 66), (66, 67), (66, 70), (66, 71), (66, 91), (67, 42), (67, 62), (67, 66), (67, 67), (67, 68), (67, 72), (67, 92), (68, 43), (68, 63), (68, 67), (68, 68), (68, 69), (68, 73), (68, 93), (69, 44), (69, 64), (69, 68), (69, 69), (69, 70), (69, 74), (69, 94), (70, 45), (70, 65), (70, 66), (70, 69), (70, 70), (70, 75), (70, 95), (71, 46), (71, 51), (71, 66), (71, 71), (71, 72), (71, 75), (71, 96), (72, 47), (72, 52), (72, 67), (72, 71), (72, 72), (72, 73), (72, 97), (73, 48), (73, 53), (73, 68), (73, 72), (73, 73), (73, 74), (73, 98), (74, 49), (74, 54), (74, 69), (74, 73), (74, 74), (74, 75), (74, 99), (75, 50), (75, 55), (75, 70), (75, 71), (75, 74), (75, 75), (75, 100), (76, 51), (76, 76), (76, 77), (76, 80), (76, 81), (76, 96), (76, 101), (77, 52), (77, 76), (77, 77), (77, 78), (77, 82), (77, 97), (77, 102), (78, 53), (78, 77), (78, 78), (78, 79), (78, 83), (78, 98), (78, 103), (79, 54), (79, 78), (79, 79), (79, 80), (79, 84), (79, 99), (79, 104), (80, 55), (80, 76), (80, 79), (80, 80), (80, 85), (80, 100), (80, 105), (81, 56), (81, 76), (81, 81), (81, 82), (81, 85), (81, 86), (81, 106), (82, 57), (82, 77), (82, 81), (82, 82), (82, 83), (82, 87), (82, 107), (83, 58), (83, 78), (83, 82), (83, 83), (83, 84), (83, 88), (83, 108), (84, 59), (84, 79), (84, 83), (84, 84), (84, 85), (84, 89), (84, 109), (85, 60), (85, 80), (85, 81), (85, 84), (85, 85), (85, 90), (85, 110), (86, 61), (86, 81), (86, 86), (86, 87), (86, 90), (86, 91), (86, 111), (87, 62), (87, 82), (87, 86), (87, 87), (87, 88), (87, 92), (87, 112), (88, 63), (88, 83), (88, 87), (88, 88), (88, 89), (88, 93), (88, 113), (89, 64), (89, 84), (89, 88), (89, 89), (89, 90), (89, 94), (89, 114), (90, 65), (90, 85), (90, 86), (90, 89), (90, 90), (90, 95), (90, 115), (91, 66), (91, 86), (91, 91), (91, 92), (91, 95), (91, 96), (91, 116), (92, 67), (92, 87), (92, 91), (92, 92), (92, 93), (92, 97), (92, 117), (93, 68), (93, 88), (93, 92), (93, 93), (93, 94), (93, 98), (93, 118), (94, 69), (94, 89), (94, 93), (94, 94), (94, 95), (94, 99), (94, 119), (95, 70), (95, 90), (95, 91), (95, 94), (95, 95), (95, 100), (95, 120), (96, 71), (96, 76), (96, 91), (96, 96), (96, 97), (96, 100), (96, 121), (97, 72), (97, 77), (97, 92), (97, 96), (97, 97), (97, 98), (97, 122), (98, 73), (98, 78), (98, 93), (98, 97), (98, 98), (98, 99), (98, 123), (99, 74), (99, 79), (99, 94), (99, 98), (99, 99), (99, 100), (99, 124), (100, 75), (100, 80), (100, 95), (100, 96), (100, 99), (100, 100), (100, 125), (101, 1), (101, 76), (101, 101), (101, 102), (101, 105), (101, 106), (101, 121), (102, 2), (102, 77), (102, 101), (102, 102), (102, 103), (102, 107), (102, 122), (103, 3), (103, 78), (103, 102), (103, 103), (103, 104), (103, 108), (103, 123), (104, 4), (104, 79), (104, 103), (104, 104), (104, 105), (104, 109), (104, 124), (105, 5), (105, 80), (105, 101), (105, 104), (105, 105), (105, 110), (105, 125), (106, 6), (106, 81), (106, 101), (106, 106), (106, 107), (106, 110), (106, 111), (107, 7), (107, 82), (107, 102), (107, 106), (107, 107), (107, 108), (107, 112), (108, 8), (108, 83), (108, 103), (108, 107), (108, 108), (108, 109), (108, 113), (109, 9), (109, 84), (109, 104), (109, 108), (109, 109), (109, 110), (109, 114), (110, 10), (110, 85), (110, 105), (110, 106), (110, 109), (110, 110), (110, 115), (111, 11), (111, 86), (111, 106), (111, 111), (111, 112), (111, 115), (111, 116), (112, 12), (112, 87), (112, 107), (112, 111), (112, 112), (112, 113), (112, 117), (113, 13), (113, 88), (113, 108), (113, 112), (113, 113), (113, 114), (113, 118), (114, 14), (114, 89), (114, 109), (114, 113), (114, 114), (114, 115), (114, 119), (115, 15), (115, 90), (115, 110), (115, 111), (115, 114), (115, 115), (115, 120), (116, 16), (116, 91), (116, 111), (116, 116), (116, 117), (116, 120), (116, 121), (117, 17), (117, 92), (117, 112), (117, 116), (117, 117), (117, 118), (117, 122), (118, 18), (118, 93), (118, 113), (118, 117), (118, 118), (118, 119), (118, 123), (119, 19), (119, 94), (119, 114), (119, 118), (119, 119), (119, 120), (119, 124), (120, 20), (120, 95), (120, 115), (120, 116), (120, 119), (120, 120), (120, 125), (121, 21), (121, 96), (121, 101), (121, 116), (121, 121), (121, 122), (121, 125), (122, 22), (122, 97), (122, 102), (122, 117), (122, 121), (122, 122), (122, 123), (123, 23), (123, 98), (123, 103), (123, 118), (123, 122), (123, 123), (123, 124), (124, 24), (124, 99), (124, 104), (124, 119), (124, 123), (124, 124), (124, 125), (125, 25), (125, 100), (125, 105), (125, 120), (125, 121), (125, 124), (125, 125)]

E23DT = []
for e in E3DT 
    if !(e[1] == e[2])
        push!(E23DT, [e[1], e[2]])
    end
end

G3DTorus = digraph(V3DT,E23DT)

An = homology.A(G3DTorus,3)
#@time homology.A(G3DTorus,3)
#homology.restrictedMatrixV2(An,3)
On = homology.O_n(An,3)

ODiff = homology.O_diffV3(On,3)

#println(ODiff[3])
@time homology.snfDiagonal(ODiff[3])
#@btime homology.pathHomologyV2(G3DTorus,1)
#sfH3DTorus  = homology.pathHomologyV2(G3DTorus,3)
#println("3D Torus: ", sfH3DTorus)
#homology.printSNF(sfH3DTorus)
#@btime homology.pathHomologyV2(G3DTorus,3)
# Example Moibus Strip 
VM = [1,2,3,4,5,6]
EM = [[1,2], [2,3], [3,1], [2,4], [4,3], [3,5], [5,4], [6,4], [5,6], [6,5]]
GM = digraph(VM,EM)
#homology.pathHomologyV2(GM,5)
#   2.132485 seconds (44.81 M allocations: 6.147 GiB, 10.94% gc time)

VH = [1,2,22,3,33,4]
EH = [[1,2], [2,3], [3,4], [1,22], [22,33],[33,4]]
hexagonUP = digraph(VH,EH)
#println(homology.pathHomologyV2(hexagonUP,5))

ETH = [[1,2], [2,3], [3,4], [1,22], [22,33],[33,4],[1,3],[2,4],[1,33],[22,4],[1,4]]
hexagonUPT = digraph(VH,ETH)
#println(homology.pathHomologyV2(hexagonUPT,5))

#V = [1,2,3]
#E = [[1,2],[2,3],[3,1]]
#E2 = [[1,2],[3,2],[3,1]]
#g1 = digraph(V,E)
#g2 = digraph(V,E2)
#println(homology.pathHomologyV2(g1,5))
#println(homology.pathHomologyV2(g2,5))

V = [1,2,3,4,5,6]
E = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 5], [4, 5], [4, 6], [5, 6]]
V2 = [1,2,3]
E2 = [[1,2],[1,3],[1,2]]
g2 = digraph(V2,E2)
#println(homology.pathHomologyV2(g2,5))

V = [1,2,3,4,5]
E = [[1,2],[2,3],[3,1],[1,4],[2,4],[3,4],[1,5],[3,5],[2,5]]
G = digraph(V,E)
H = [V,E]
#println("Hypergraph: ", homology.HyperPathHomologyV2(H,2,5))
#H = [[1,2,3,4,5],[[2,1],[4,1],[3,1],[2,3],[4,2],[3,4],[5,2],[4,5],[3,5]]]
#G = digraph(H[1],H[2])
#println("Suspension: ", homology.pathHomologyV2(G,5))
#println(homology.A(G,4))

#=====================================================================================================================
HYPERGRAPH EXAMPLES 
======================================================================================================================#
V = [1,2,3]
E = [[1,2,3], [1,2], [2,3]]
E2 = [[1,2],[2,3]]

V3 = [1,2,3,4]
E3 = [[1,2,3],[1,2,4],[1,3,4],[2,3,4]]

V4 = [1,2,3,4,5]
E4 = [[1,2,3,4,5], [1,2,3],[1,2,4],[1,3,4],[2,3,4]]

V5 = [1,2,3,4,5,6,7]
E5 = [[1,2,3,4,5],[1,2,3,4,6,7],[1,2,3],[1,3,4],[1,2,4]]

H = [V,E]
#println(homology.HyperPathHomology(H,2,5))

H2 = [V,E2]
#println(homology.HyperPathHomology(H2,2,5))

H3 = [V3,E3]
#println(homology.HyperPathHomology(H3,2,5))
 
H4 = [V4,E4] 
#println(homology.HyperPathHomology(H4,4,5))

H5 = [V5,E5]
#println(homology.HyperPathHomology(H4,2,5))

# function to remove edges of degree n 

function cleanDegn(V,E)
    for v in V
        incoming = false
        outgoingEdgeSet = [] 
        outgoingVertexSet = [] 
        for e in E 
            if v == e[2] # check if the vertex has an incoming edge
                incoming = true
                break
            else
                if v == e[1] 
                    push!(outgoingEdgeSet, e)
                    push!(outgoingVertexSet, e[2])
                end
            end
        end
        outgoingCond = false
        if incoming == false
            n = length(outgoingVertexSet)
            if n > 1
                for i in 1:n
                    outgoingVertex = outgoingVertexSet[i]
                    tempSet = copy(outgoingVertexSet)
                    for edge in E
                        if outgoingVertex == edge[1]
                            if edge[2] in tempSet
                                deleteat!(tempSet, findall(x->x==edge[2],tempSet))

                            end
                        end
                    end
                    deleteat!(tempSet, findall(x->x==outgoingVertex,tempSet))
                    if isempty(tempSet) == true
                        outgoingCond = true
                        break
                    end
                end

                if outgoingCond == true
                    deleteat!(V, findall(x->x==v,V))
                    for e in outgoingEdgeSet
                        deleteat!(E, findall(x->x==e,E))
                    end
                end
            end
        end
    end
    return V,E
end

# function to determine if two vertices are an edge
function isEdge(e, E)
    for edge in E
        if (e[1] == edge[1]) && (e[2] == edge[2])
            return true
        end

        if (e[1] == edge[2]) && (e[2] == edge[1])
            return true
        end
    end
    return false
end

# function to determine if an two vertices are a semi-edge
function isSemiEdge(e, E)
    c = e[1]
    cList = []
    b = e[2]
    bList = []

    for edge in E
        if c == edge[1]
            push!(cList,edge)
        end
        if b == edge[2]
            push!(bList, edge)
        end
    end

    for cEdge in cList
        for bEdge in bList
            if cEdge[2] == bEdge[1]
                return true
            end
        end
    end
    return false
end

# function to partition vertices into the graphs connected componets
function connectedComponets(V,E)
    partition = []
    for v in V
        vList = [] 
        push!(vList, v)
        for e in E
            if v == e[1]
                push!(vList, e[2])
            end
            if v == e[2]
                push!(vList, e[1])
            end
        end
        if isempty(partition) == true
            push!(partition, vList)
        else
            for list in partition
                intersection = intersect(list, vList)
   
                if isempty(intersection) == false
                    pos = findall(x->x==list, partition)
                    combine = union(list, vList)
                    comb = sort(combine)
                    if comb in partition 
                    else
                        replace(partition, list =>comb)
                    end
                    #partition[pos] = combine
                else
                    if sort(vList) in partition
                    else
                        push!(partition, vList)
                    end
                end
            end
        end
    end
    return partition
end

# function to determine if two vertices are in the same componet
function sameComponent(e, E)
    b = e[1]
    c = e[2]
    for edgeSet in E
        if (b in edgeSet) && (c in edgeSet)
            return true
        end
    end
    return false
end

# function to remove edges of degree n 
function cleanDegn2(V,E)
    for v in V
        incoming = false
        incomingEdgeSet = [] 
        incomingVertexSet = [] 
        for e in E 
            if v == e[1] # check if the vertex has an outgoing edge
                incoming = true
                break
            else
                if v == e[2] 
                    push!(incomingEdgeSet, e)
                    push!(incomingVertexSet, e[1])
                end
            end
        end
        incomingCond = false
        if incoming == false
            n = length(incomingVertexSet)
            if n > 1
                for i in 1:n
                    incomingVertex = incomingVertexSet[i]
                    tempSet = copy(incomingVertexSet)
                    for edge in E
                        if incomingVertex == edge[2]
                            if edge[1] in tempSet
                                deleteat!(tempSet, findall(x->x==edge[1],tempSet))
                            end
                        end
                    end
                    deleteat!(tempSet, findall(x->x==incomingVertex,tempSet))
                    if isempty(tempSet) == true
                        incomingCond = true
                        break
                    end
                end

                if incomingCond == true
                    deleteat!(V, findall(x->x==v,V))
                    for e in incomingEdgeSet
                        deleteat!(E, findall(x->x==e,E))
                    end
                end
            end
        end
    end
    return V,E
end

function cleanEdges(V,E)
    for e in E 
        branch = true

        incoming = false
        incomingCount = 0
        outgoing = false
        outgoingCount = 0
        for e2 in E
            if !(e == e2)
                if (e[2] == e2[1])
                    outgoing = true
                    outgoingCount = outgoingCount + 1
                elseif (e[1] == e2[2])
                    incoming = true
                    incomingCount = incomingCount + 1
                elseif (e[1] == e2[1])
                    incoming = true
                    incomingCount = incomingCount + 1
                elseif (e[2] == e2[2])
                    outgoing = true
                    outgoingCount = outgoingCount + 1
                end
            end
        end
        if (outgoing == true) && (incoming == true)
            branch = false
        elseif (outgoing == true) && (incoming == false)
        elseif (outgoing == false) && (incoming == true)
            if incomingCount > 1 
                branch = false
            end
        end

        if branch == true
            deleteat!(E, findall(x->x==e,E))
        end
    end

    for v in V
        remove = true
        for e in E
            if (v == e[1])
                remove = false
                break
            elseif (v == e[2])
                remove = false
                break 
            end
        end
        if remove == true
            deleteat!(V, findall(x->x==v,V))
        end
    end
    return V,E
end

function cleanDeg1plus1(V,E)
    changeHomology1 = []
    changeHomology0 = []
    #VertexPartition = connectedComponets(V,E)
    for v in V
        incomingEdges = []
        incomingVertex = [] 
        outgoingEdges = []
        outgoingVertex = []
        for e in E
            if v == e[1]
                push!(outgoingEdges,e)
                push!(outgoingVertex,e[2])


                if length(outgoingEdges) > 1
                    break
                end
            end
            if v == e[2]
                push!(incomingEdges, e)
                push!(incomingVertex,e[1])

                if length(incomingEdges) > 1
                    break
                end
            end
        end

        if isempty(incomingVertex) == true
            break
        end
        if isempty(outgoingVertex) == true
            break
        end

        if !(incomingVertex[1] == outgoingVertex[1])

            deleteat!(E, findall(x->x==incomingEdges[1],E))
            deleteat!(E, findall(x->x==outgoingEdges[1],E))

            deleteat!(V, findall(x->x==v,V))

            semiE = isSemiEdge([incomingVertex[1],outgoingVertex[1]], E)
            edge = isEdge([incomingVertex[1],outgoingVertex[1]], E)

            if (semiE == true) || (edge == true)
                push!(changeHomology0, false)
                push!(changeHomology1, false)
            end

            if (semiE == false) && (edge == false)
                VertexPartition = connectedComponets(V,E)
                connected = sameComponent([incomingVertex[1],outgoingVertex[1]],VertexPartition)

                if connected == true
                    push!(changeHomology0, false)
                    push!(changeHomology1, true)
                else
                    push!(changeHomology0, true)
                    push!(changeHomology1, false)
                end
            end
        end
    end
    return V, E
end