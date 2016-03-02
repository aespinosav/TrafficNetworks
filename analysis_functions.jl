#This file contains functions for getting the total cost and other
#useful things for analysing simulation data.

"""
Returns the total cost of the network for every demand
value in the demand range.

TA problem has to have been solved
"""
function total_cost(rn::RoadNetwork)
    m = num_edges(rn)
    num_sols = size(rn.flows)[2]
    num_qs = length(rn.demand_range)
    costs = [dot(rn.a,rn.flows[:,i]) + dot(rn.flows[:,i],diagm(rn.b)*rn.flows[:,i]) for i in 1:num_qs]
end

"""
Returns the actual cost of each edge (cost * flow) in the same format as the
flows (rows correspond to edges, columns to demand level)
"""
function edge_costs(rn::RoadNetwork)
    m = num_edges(rn)
    num_sols = size(rn.flows)[2]
    num_qs = length(rn.demand_range)
    
    costs = rn.flows .* (rn.a .+ rn.b .* rn.flows)
end