using JuMP, Gurobi

function make_demand_vector(rn, od, demand)
    n = num_nodes(rn.g)
    d = spzeros(n)
    
    indices = find(od)
    for k in indices
        i, j = ind2sub(od, k)
        if od[i,j] > 0
            d[i] = -od[i,j]*demand
            d[j] =  od[i,j]*demand
        end
    end
    d
end


function stap_solve(rn::RoadNetwork, od, demand_range::Array{Float64,1}; regime="UE")
    m = num_edges(rn.g)
    n = num_nodes(rn.g)
    a = rn.a
    b = rn.b
    inc_mat = incidence_matrix(rn.g)
    d = make_demand_vector(rn, od, demand_range[1])
    
    STAP = Model(solver=Gurobi.GurobiSolver())
    
    @variable(STAP, x[1:m])
    @constraint(STAP, positive_flow_const, x.>=0)
    @constraint(STAP, conservation_const[i=1:n], AffExpr(x, inc_mat[i,:], 0.0) == d[i])
    
    beckmann = QuadExpr(x, x, 0.5*b, AffExpr(x, a, 0.0))
    sys_cost = QuadExpr(x, x, b, AffExpr(x, a, 0.0))
    
    if regime=="UE"
        @objective(STAP, Min, beckmann)
    elseif regime=="SO"
        @objective(STAP, Min, sys_cost)
    else
        error("Regime not recognised... try 'SO' or 'UE'\n")
    end
        
    show(STAP)
    
    sols = Array{Float64,1}[]
    solve(STAP) 
    push!(sols, getvalue(x))
    for q in demand_range[2:end]
        d = make_demand_vector(rn, od, q)
        for i in 1:n
            JuMP.setRHS(conservation_const[i], d[i])
        end
        solve(STAP)
        push!(sols, getvalue(x))
    end
    sols = hcat(sols...)
end
stap_solve(rn::RoadNetwork, od, demand::Float64; Regime="UE") = stap_solve(rn, od, [demand], regime=Regime)[:]
