# Types
##########################
##########################

"""
Road network object, has a graph, an OD matrix with the associated flows 
between the O-D pair, and a and b parameters for affine cost functions.

Can be constructed giving either a graph or and adjacency matrix as the first parameter:
 
RoadNetwork(g, a, b, OD)
RoadNetwork(A, a, b, OD)

I should probably define a convert function from matrices to graphs to make code neater...
"""
type RoadNetwork
    g::Graph
    a::Array{Float64,1}
    b::Array{Float64,1}
    OD::Array{Float64,2}
    demand_range::Array{Float64,1}
    flows_ue::Array{Float64,2}
    flows_so::Array{Float64,2}
    flows_other::Array{Float64,2}
end

function RoadNetwork(g::Graph, a::Array{Float64,1}, b::Array{Float64,1}, OD::Array{Int64,2})
    RoadNetwork(g, a, b, OD, Array{Float64,1}(), Array{Float64,2}(), Array{Float64,2}(),Array{Float64,2}())
end

RoadNetwork(A::Array{Int,2}, a::Array{Float64,1}, b::Array{Float64,1}, OD::Array{Int64,2}) =
begin
RoadNetwork(Graph(A), a, b, OD, Array{Float64,1}(), Array{Float64,2}(), Array{Float64,2}(), Array{Float64,2}())
end

function show(io::IO, rn::RoadNetwork)
    output = "RoadNetwork:\nNodes - $(num_nodes(rn.g)) \nEdges - $(num_edges(rn.g))\n    OD\n$(rn.OD)"
    print(io, output)
end

# Functions
##########################
##########################

num_nodes(rn::RoadNetwork) = num_nodes(rn.g)
num_edges(rn::RoadNetwork) = num_edges(rn.g)

"""
Replaces OD matrix in a RoadNetwork rn.
"""
function replace_OD_matrix!(rn::RoadNetwork, OD_new::Array{Float64,2})
    rn.OD = OD_new
end

#question: should I attach the flows to the RoadNetwork object??