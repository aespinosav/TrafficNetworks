    # a destination based variable choice
# Still need to test on network from ensemble.

using TrafficNetworks, JuMP, Gurobi, SkeletonCities


#function make_demands_mat(g, od, demands)
#   n = num_nodes(g)
#    num_flows = length(demands)
#    d_mat = spzeros(n, num_flows)
#    
#    counter = 0
#    for i in 1:n
#        for (j,k) in enumerate(find(od[:,i]))
#            counter += 1
#            d_mat[i,counter] = sum(demands[counter]) # no self demand...
#            d_mat[k,counter] = -demands[counter][j] #because demands is an array of arrays...
#        end
#    end
#    d_mat
#end


"""
Returns the node index of the destination associated to
the flow numbering passed as `k` 
"""
function dest_nodes_of_flows(od)
    pair_indices = find(od')
    num_pairs = length(pair_indices)

    pair_origin_idxs = zeros(Int64, num_pairs)
    pair_destination_idxs = zeros(Int64, num_pairs)
    for i in 1:num_pairs
        pair_destination_idxs[i], pair_origin_idxs[i] =  ind2sub(od', pair_indices[i])
    end
    
    dest_indices = sort(unique(pair_destination_idxs))
end

"""
Auxiliary function to make rest of code tidier
"""
function num_flows(od::SparseMatrixCSC)
    length(dest_nodes_of_flows(od))
end

function dest_tree(od::SparseMatrixCSC)
    dest_indices = dest_nodes_of_flows(od)
    tree = Dict{Int64, Array{Int64,1}}()
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
    
    #stap_multi = Model(solver=Gurobi.GurobiSolver())
    stap_multi = Model(solver=Gurobi.GurobiSolver(OutputFlag=0))

    @variable(stap_multi, x[1:m,1:nf] >= 0)
    @variable(stap_multi, rowsums[1:m] >= 0) #Aggregate link variables
    
    @constraint(stap_multi, incidence_matrix(rn.g)*x .== d_mat)
    @constraint(stap_multi, inter_var_con[i in 1:m], rowsums[i] == sum(x[i,j] for j in 1:nf) )    
    @objective(stap_multi, Min, dot(a,rowsums) + (C*rowsums'*B*rowsums)[1])
    
    solve(stap_multi)
    
    # Sanity check, but make sure they are the same dims
    #if sum(getvalue(x,2)) == getvalue(rowsums)
    #   println("Aggregates check out!") 
    #end
    getvalue(x)
end
