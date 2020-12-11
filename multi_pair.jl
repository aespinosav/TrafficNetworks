#This script does not work, might be better to just get rid of it as multi_pair_e.jl now seems to work

function make_demand_vectors(g, od, demands)
    n = num_nodes(g)
    m = num_edges(g)
    
    d_vecs = SparseVector{Float64,Int64}[]
    
    indices = find(od)
    flow_counter = 0
    for k in indices
        i, j = ind2sub(od, k)
        if od[i,j] > 0
            flow_counter += 1
            d = demands[flow_counter]
            
            d_vec = spzeros(n)
            
            d_vec[i] = -od[i,j]*d
            d_vec[j] =  od[i,j]*d
            
            push!(d_vecs, d_vec)
        end
    end
    d_vecs
end


function multipair_objective(x, k, num_flows, a, b)
    m = length(a)
    B = diagm(b)

    if k == 1
        y = sum([getvalue(x[i]) for i in 2:num_flows])
    elseif k < num_flows
        y = sum([getvalue(x[i]) for i in 1:k-1]) + sum([getvalue(x[i]) for i in k+1:num_flows])
    else
        y = sum([getvalue(x[i]) for i in 1:num_flows])
    end 
    
    obj = a'*x[k] + 0.5*x[k]'*B*x[k] + y'*B*x[k]
    obj[1]
end

function multipair_stap(rn, od, demands; tolerance=1e-6, max_iters=10)

    m = num_edges(rn.g)
    n = num_nodes(rn.g)
    inc_mat = incidence_matrix(rn.g)
    a = rn.a
    b = rn.b
    num_flows = length(demands)
    indices = find(od)
    
    x_old = [zeros(m) for k in 1:num_flows]
    STAP = [Model(solver=GurobiSolver()) for i in 1:num_flows]
    x = Array{Array{JuMP.Variable,1},1}(num_flows)
    
    d_vecs = make_demand_vectors(rn.g, od, demands)
    
    for k in 1:num_flows
        
        ### Set up STAP ###    
        
        x[k] = @variable(STAP[k], [1:m], lowerbound=0)
        @constraint(STAP[k], AffExpr(x[k], inc_mat[k,:], 0.0) .== d_vecs[k])

        ### Init values ###
        
        od_k = od_matrix_from_pair(rn.g, ind2sub(od, k)) 
        
        x_init = ta_solve(rn, od_k, demands[k])
        setvalue(x[k], x_init)
        x_old[k] = x_init
    end 
    
    iters = 0
    last_delta = 100
    x_agg_new = zeros(m)
    for k in 1:num_flows
        x_agg_new += x_old[k]
    end
    
    while (iters < max_iters) && (last_delta > tolerance)
        iters += 1
        
        for k in 1:num_flows
            @objective(STAP[k], Min, multipair_objective(x[k], k, num_flows, a, b))
            solve(STAP[k])
            x_new[k] = getvalue(x[k])
            x_agg_new += x_new[k]
        end
        
        last_delta = norm(x_agg_new - x_agg_old)
        x_agg_old = x_agg_new        
    end
         
    sols = [getvalue(x[k]) for k in 1:num_flows]
end
