#In this script we implement eddie's suggested version (a destination based variable choice)
#Still need to test on network from ensemble.

using TrafficNetworks, JuMP, Gurobi, SkeletonCities


"""
Returns the node index of the destination associated to
the flow numbering passed as `k` 
"""
function dest_nodes_of_flows(od)
    pair_indices = find(od')
    num_pairs = length(pair_indices)

    pair_origin_idxs = zeros(Int, num_pairs)
    pair_destination_idxs = zeros(Int, num_pairs)
    for i in 1:num_pairs
        pair_destination_idxs[i], pair_origin_idxs[i] =  ind2sub(od', pair_indices[i])
    end
    
    dest_indices = sort(unique(pair_destination_idxs))
end

function num_flows(od)
    length(dest_nodes_of_flows(od))
end

function dest_tree(od)
    dest_indices = dest_nodes_of_flows(od)
    tree = Dict{Int, Array{Int,1}}()
    for (i, k) in enumerate(dest_indices)
        # i := num flows
        # k := destination node index
        origin_list = find(od[:,k])
        tree[k] = origin_list
    end
    tree
end

function origins_per_destination(tree)
    [length(tree[k]) for k in dest_nodes_of_flows]
end



"""
Retruns a demand matrix where each column has, for each flow (labeled by detination, in ascending index order),
the source and sink values needed for the flow conservation constraints.

    `make_demands_mat(g, od, demands)`
"""
function make_demands_mat(g, od, demands)    
    
    n = num_nodes(g)
    
    nf = num_flows(od)   
    dest_indices = dest_nodes_of_flows(od)
    flow_tree = dest_tree(od)

    d_mat = spzeros(n, nf)    
    for (i, k) in enumerate(dest_indices)
        origins = flow_tree[k]
        d_mat[k,i] = sum(demands[i])#destination condition
        for (j,o) in enumerate(origins)
            d_mat[o,i] = -demands[i][j]#origin condition
        end
    end
    
    d_mat
end


function multi_pair_stap(rn, od, demands; regime="UE")
    C = (regime == "UE" ? 0.5 : 1.0)
    
    n = num_nodes(rn)
    m = num_edges(rn)
    
    B = diagm(rn.b)
    a = rn.a
    
    nf = num_flows(od)
    flow_tree = dest_tree(od)
    
    d_mat = make_demands_mat(rn.g, od, demands)
    
    stap_multi = Model(solver=Gurobi.GurobiSolver())

    @variable(stap_multi, x[1:m,1:nf] >= 0)
    @variable(stap_multi, rowsums[1:m] >= 0) #Aggregate link variables
    
    
    @constraint(stap_multi, incidence_matrix(rn.g)*x .== d_mat)
    @constraint(stap_multi, sum(x, 2)  .<= 0) #capacity constraint
    @constraint(stap_multi, inter_var_con[i in 1:m], rowsums[i] == sum(x[i,j] for j in 1:nf) )    
    @objective(stap_multi, Min, dot(a,rowsums) + (C*rowsums'*B*rowsums)[1])
    
    solve(stap_multi)
    
    #sanity check, but make sure they are the same dims
    #if sum(getvalue(x,2)) == getvalue(rowsums)
    #   println("Aggregates check out!") 
    #end
    getvalue(x)
end

