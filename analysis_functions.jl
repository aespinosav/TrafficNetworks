#This file contains functions for getting the total cost and other
#useful things for analysing simulation data.

"""
Returns the total cost of the network for every demand
value in the demand range.

TA problem has to have been solved
"""
function total_cost(rn::RoadNetwork; regime="UE")
    if regime == "UE"
        flows = rn.flows_ue
    else
        flows = rn.flows_so
    end
    
    m = num_edges(rn)
    num_sols = size(flows)[2]
    num_qs = length(rn.demand_range)
    costs = [dot(rn.a,flows[:,i]) + dot(flows[:,i],diagm(rn.b)*flows[:,i]) for i in 1:num_qs]
end

"""
Returns the actual cost of each edge (cost * flow) in the same format as the
flows (rows correspond to edges, columns to demand level)
"""
function edge_costs(rn::RoadNetwork; regime="UE")
    m = num_edges(rn)
    
    if regime=="UE"
        flows = rn.flows_ue
    else
        flows = rn.flows_so
    end
    
    num_sols = size(flows)[2]
    num_qs = length(rn.demand_range)
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

    poa = rn.flows_ue ./ rn.flows_so
end