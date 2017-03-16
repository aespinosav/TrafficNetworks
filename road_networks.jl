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

"""
Constructs a RoadNetwork from graph 'g' with param vectors 'a' and 'b' (lin cost functions).
and an 'od_pair' (given as) either a tuple or an array
"""
function RoadNetwork(g::Graph, a::Array{Float64,1}, b::Array{Float64,1}, od_pair)
    OD = od_matrix_from_pair(g, od_pair)
    RoadNetwork(g, a, b, OD)
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
incidence_matrix(rn::RoadNetwork) = incidence_matrix(rn.g)

"""
Replaces OD matrix in a RoadNetwork rn.
"""
function replace_OD_matrix!(rn::RoadNetwork, OD_new::Array{Float64,2})
    rn.OD = OD_new
end

"""
Generates OD matrix for graph 'g' for a single OD pair given as a 2-element array 'od_pair'

This should also be changed to sparse matrices once the rest of the code gets a rehaul.
"""
function od_matrix_from_pair(g, od_pair)
    n = num_nodes(g)
    OD = zeros(Int, n, n)
    OD[od_pair[1], od_pair[2]] = 1
    OD
end

#question: should I attach the flows to the RoadNetwork object??
#
# I think I should define a new "Assignment" object since it is getting a bit clunky to 
# keep putting everything into the RoadNetwork. That way the ta_solve functions can take a RoadNetwork 
# onject and return an assignement, that focuses on the flows and costs and all the data of the assignment
# whilst keeping the RoadNetwork objec as minimal as needed. This would meand that the OD matrix just gives 
# pairs of origins and destinations and giving these demands ranges would just have to be through the number
# of distinct flows.
#
# Doing this also means that once the traffic assignment has been solved, we deal only with data from the simulation
# and not with objects that can be changed by the solver functions. Saving the assignment would therefore be quite
# simple, we have the header for the details  of the simulation and a table for the flows. Which can then be imported
# into any  language easily for processing etc, in a simple way.
#
# All this should be done in a new branch!
