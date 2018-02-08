#This file contains functions for getting the total cost and other
#useful things for analysing simulation data.

#"""
#Returns the total cost of the network for every demand
#value in the demand range.
#
#TA problem has to have been solved
#"""
#function total_cost(rn::RoadNetwork, flows, demand_range)
#    m = num_edges(rn)
#    num_sols = size(flows)[2]
#    num_qs = length(demand_range)
#    costs = [dot(rn.a,flows[:,i]) + dot(flows[:,i],diagm(rn.b)*flows[:,i]) for i in 1:num_qs]
#end

"""
Returns the total cost (scalar value) of the assignment given by flow_vector
"""
function total_cost(rn::RoadNetwork, flow_vector::Array{Float64,1})
    total_cost = rn.a'*flow_vector + flow_vector'*diag(rn.b)*flow_vector
end

"""
Returns a 1D array of the total network cost for each corresponding flow vector (given as the columns of flow_vector_array)
"""
function total_cost(rn:RoadNetwork, flow_vector_array::Array{Float64,2})
    cols = size(flow_vector_array)[2]
    total_costs_array = Array{Float64}(cols)
    for i in 1:cols
        total_costs_array[i] = total_cost(rn, flow_vector_array[:,i][:])
    end
    total_costs_array
end

"""
Returns a 1D cost vector of the costs of each for the flow assignment given by flow vector
"""
function costs(rn::RoadNetwork, flow_vector::Array{Float64,1})
    cost_vector = rn.a + (rn.b .* flow_vector)
end

"""
Returns a 2D array where each column is the cost vector for the corresponding flow assignment
given by the column of the 2D array flow_vector_array
"""
function costs(rn::RoadNetwork, flow_vector_array::Array{Float64,2})
    cols = size(flow_vector_array)[2]
    cost_vector_array = Array{Float64}(num_edges(rn.g), cols)
    for i in 1:cols
        cost_vector_array[:,i] = costs(rn, flow_vector_array[:,i][:])
    end 
    cost_vector_array
end




#
# I dont think I use this terminology anymore, plus same as function that has been commented out below...
# 
# RoadNetwork object has changed since this function was last used
#"""
#Returns the marginal cost of each edge (cost * flow) in the same format as the
#flows (rows correspond to edges, columns to demand level)
#"""
#function marginal_edge_costs(rn::RoadNetwork; regime="UE")
#    m = num_edges(rn)
#
#    if regime=="UE"
#        flows = rn.flows_ue
#    else
#        flows = rn.flows_so
#    end
#
#    num_sols = size(flows)[2]
#    num_qs = length(rn.demand_range)
#    costs = (rn.a .+ rn.b .* flows)
#end

#function tot_cost(rn::RoadNetwork, flows)
#    costs = sum(rn.a .+ rn.b .* flows)
#end


#
# The following function is outdated since RoadNetwork no longer stores the assignmnets 
# therefore I have commented out the function (should probably be removed)
#
#"""
#Returns the actual cost of each edge (cost * flow) in the same format as the
#flows (rows correspond to edges, columns to demand level)
#"""
#function edge_costs(rn::RoadNetwork; regime="UE")
#    m = num_edges(rn)
#    
#    if regime=="UE"
#        flows = rn.flows_ue
#    else
#        flows = rn.flows_so
#    end
#    
#    num_sols = size(flows)[2]
#    num_qs = length(rn.demand_range)
#    costs = flows .* (rn.a .+ rn.b .* flows)
#end
#
#"""
#Returns the price of anarchy for the demand range that the TA on the network has been solved for.
#The function checks that the network has been solved for both UE and SO.
#
#usage:
#price_of_anarchy(rn)
#"""
#function price_of_anarchy(rn::RoadNetwork)
#    #We check that flows for both cases exist
#    if size(rn.flows_so) != size(rn.flows_ue)
#        error("Mismatch in sizes for UE and SO flows")
#    end
#    poa = rn.flows_ue ./ rn.flows_so
#end
