#This file contains functions for getting the total cost and other
#useful things for analysing simulation data.

"""
Returns the total cost of the network for every demand
value in the demand range.

TA problem has to have been solved
"""
function total_cost(ta::TrafficAssignment)
    m = num_edges(ta.rn)
    #num_sols = size(flows)[2]
    num_qs = length(ta.demand_range)

    costs = [dot(ta.rn.a, Array(ta.flows[i, 2:end])'[:]) + dot(Array(ta.flows[i, 2:end])'[:],diagm(ta.rn.b)*Array(ta.flows[i, 2:end])'[:]) for i in 1:num_qs]
end

"""
Returns the marginal cost of each edge (cost * flow) in the same format as the
flows (rows correspond to edges, columns to demand level)
"""
function marginal_edge_costs(ta::TrafficAssignment)
    flows = Array(ta.flows[2:end])'
    costs = (ta.rn.a .+ (ta.rn.b .* flows))
end


"""
Returns the actual cost of each edge (cost * flow) in the same format as the
flows (rows correspond to edges, columns to demand level)
"""
function edge_costs(ta::TrafficAssignment)
    m = num_edges(rn)
    
    flows = Array(ta.flows[2:end])'
    costs = flows .* (rn.a .+ rn.b .* flows)
end

"""
Returns the price of anarchy for the demand range that the TA on the network has been solved for.
The function checks that the network has been solved for both UE and SO.

usage:
price_of_anarchy(rn)
"""
function price_of_anarchy(rn::RoadNetwork)
    #We check that flows for both cases exist
    if size(rn.flows_so) != size(rn.flows_ue)
        error("Mismatch in sizes for UE and SO flows")
    end
    poa = rn.flows_ue ./ rn.flows_so
end