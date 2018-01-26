#In this file we implement the RoadNetwork object and its pertaining functions.
#It has changed a bit since the beginning and it is not clear if this type is actually just cluttering up this library.

# Types
"""
Road network object, has a graph, an OD matrix with the associated flows 
between the O-D pair, and a and b parameters for affine cost functions.

Can be constructed giving either a graph or and adjacency matrix as the first parameter:
 
RoadNetwork(g, a, b)
RoadNetwork(A, a, b)
"""
type RoadNetwork
    g::Graph
    a::Array{Float64,1}
    b::Array{Float64,1}
end

RoadNetwork(A::Array{Int,2}, a::Array{Float64,1}, b::Array{Float64,1}) =
begin
    RoadNetwork(Graph(A), a, b)
end

RoadNetwork(A::AbstractSparseMatrix, a::Array{Float64,1}, b::Array{Float64,1}) =
begin
    RoadNetwork(Graph(A), a, b)
end

function show(io::IO, rn::RoadNetwork)
    output = "RoadNetwork:\nNodes - $(num_nodes(rn.g)) \nEdges - $(num_edges(rn.g))\n"
    print(io, output)
end

#Functions for RoadNetwork object
num_nodes(rn::RoadNetwork) = num_nodes(rn.g)
num_edges(rn::RoadNetwork) = num_edges(rn.g)
incidence_matrix(rn::RoadNetwork) = incidence_matrix(rn.g)

"""
Returns an array of N OD pairs for a given graph g (uniformly chosen at random)
"""
function od_pairs(g::Graph, N)

    od_pair_array = Any[]
    for i in 1:N

        origin = rand(1:num_nodes(g))
        destination = rand(1:num_nodes(g))
        while destination == origin
            destination = rand(1:num_nodes(g))
        end

        push!(od_pair_array, (origin, destination))
    end
    od_pair_array
end

"""
Returns an array of N OD pairs for a given road network rn  (uniformly chosen at random).

Same as function for just graphs but takes road networks and extracts the graph object
"""
od_pairs(rn::RoadNetwork, N=1) = od_pairs(rn.g, N)

"""
Generates OD matrix for graph 'g' for a single OD pair given as a 2-element tuple 'od_pair'.
Returns a sparse matrix.
"""
function od_matrix_from_pair(g::Graph, od_pair::Tuple{Int64,Int64})
    n = num_nodes(g)
    OD = spzeros(Int64, n, n)
    OD[od_pair[1], od_pair[2]] = 1
    OD
end

"""
If given a list of od pairs as an array of tuples 'od_pairs', will construct
sparse OD matrix.
"""
function od_matrix_from_pair(g::Graph, od_pairs::Array{Tuple{Int64,Int64},1})
    n = num_nodes(g)
    OD = spzeros(Int64, n, n)
    for od in od_pairs
        OD[od[1], od[2]] = 1
    end
    OD
end


"""
Same as other od_matrix_from_pair function but takes a RoadNetwork instead of a Graph.
"""
function od_matrix_from_pair(rn::RoadNetwork, od_pair::Tuple{Int64,Int64})
    od_matrix_from_pair(rn.g, od_pair)
end
"""
Same as other function that makes an OD matrix for a netwrk from an array of od pairs.
This function takes a RoadNetwork instead of a Graph.
"""
function od_matrix_from_pair(rn::RoadNetwork, od_pairs::Array{Tuple{Int64,Int64},1} )
    od_matrix_from_pair(rn.g, od_pairs)
end

"""
Same as od_matrix_from_pair but returns a full matrix.
"""
function od_matrix_from_pair_non_sparse(g, od_pair)
    OD = od_matrix_from_pairs(g, od_pair)
    full(OD)
end

function node_positions(rn::RoadNetwork)
    node_positions(rn.g)
end
