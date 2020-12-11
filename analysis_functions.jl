#This file contains functions for getting the total cost and other
#useful things for analysing simulation data.


"""
Returns the total cost (scalar value) of the assignment given by flow_vector

    total_cost(rn::RoadNetwork, flow_vector::Array{Float64,1})
    
or
    
    `flow_vector_array::Array{Float64,2}`
    

NOTE: The latter method is to be used for solving TA for a range of demands. Where each demand is a column
of the flows matrix. Returns a 1D array with one element corresponding to each demand column.
"""
function total_cost(rn::RoadNetwork, flow_vector::Array{Float64,1})
    total_cost = rn.a'*flow_vector + flow_vector'*diagm(rn.b)*flow_vector
    total_cost[1]
end

"""
Returns a 1D array of the total network cost for each corresponding flow vector (given as the columns of flow_vector_array)

NOTE: This is to be used for solving TA for a range of demands. Where each demand is a column
of the flows matrix.
"""
function total_cost(rn::RoadNetwork, flow_vector_array::Array{Float64,2})
    cols = size(flow_vector_array)[2]
    total_costs_array = Array{Float64}(cols)
    for i in 1:cols
        total_costs_array[i] = total_cost(rn, flow_vector_array[:,i])
    end
    total_costs_array
end

"""
Price of Anarchy of UE assignment on RoadNetwork rn.

    `poa(rn, ue_flow, so_flow)`
"""
function poa(rn, ue_flow, so_flow)
    ue_cost, so_cost = total_cost(rn, ue_flow), total_cost(rn, so_flow)
    PoA = ue_cost / so_cost
end

"""
Returns a 1D cost vector of the costs of each for the flow assignment given by flow vector
"""
function costs(rn::RoadNetwork, flow_vector::Array{Float64,1})
    cost_vector = rn.a + (rn.b .* flow_vector)
end

"""
Returns a 2D array where each column is the cost vector for the corresponding flow assignment
given by the column of the 2D array flow_vector_array.

NOTE: This is to be used for solving TA for a range of demands. Where each demand is a column
of the flows matrix.
"""
function costs(rn::RoadNetwork, flow_vector_array::Array{Float64,2})
    cols = size(flow_vector_array)[2]
    cost_vector_array = Array{Float64}(num_edges(rn.g), cols)
    for i in 1:cols
        cost_vector_array[:,i] = costs(rn, flow_vector_array[:,i][:])
    end 
    cost_vector_array
end
