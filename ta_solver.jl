# Functions and calls to Convex.jl for solving Traffic Assignment

"""
Generates convex optimisation problem from the graph g, the Origin-Destination matrix OD,
and the demand vector q. If only one OD pair, then q is a scalar. Regime is either "UE" (default)
which is "user equilibrium" or "SO": system optimal.

This should be a general method for multiple OD pairs, but it only works for 1 OD pair for now

This function uses Convex.jl
"""
function make_ta_problem(rn::RoadNetwork, q, regime)
    
    n = num_nodes(rn)
    m = num_edges(rn)
    a = rn.a
    b = rn.b

    #num_flows = sum(map(sign, rn.OD))

    M = incidence_matrix(rn.g)

    x = Variable(m)
    d = zeros(m)
    
    flow_counter = 1
    for i in 1:n
        for j in 1:n
            if rn.OD[i,j] > 0.0
                d[i] = -q[flow_counter]*rn.OD[i,j]
                d[j] = q[flow_counter]*rn.OD[i,j]
                flow_counter += 1
            end
        end
    end

    if regime == "UE"
        cost_function = dot(a, x) + 0.5*quadform(x,diagm(b))
    elseif regime == "SO"
        cost_function = dot(a,x) + quadform(x,diagm(b))
    else
        error("Regime must be either 'UE' or 'SO'...")
    end

    # Constraints to be passes to Convex.jl problem
    eq_constraints = incidence_matrix(rn.g)*x == d
    ineq_constraints = x >= 0

    problem = minimize(cost_function, eq_constraints, ineq_constraints)

    return problem, x
end


"""
Returns solutions to the traffic assignment problem for a given range of demands: q_range.
Calls function make_ta_problem
"""
function ta_solve(rn::RoadNetwork, q_range::Array{Float64,1}, regime="UE")
    println("Will solve $regime, TA problem  for $(length(q_range)) values of demand...\n")
    
    #redirect output of Convex solver to a log file to avoid screen clutter
    originalSTDOUT = STDOUT
    f = open("log_ta_solve.txt", "w")
    redirect_stdout(f)

    sols = Array{Float64}[]
    for q in q_range
        problem, x = make_ta_problem(rn, q, regime)
        solve!(problem)
        #println("\n\nq = $q\nStatus: $(problem.status)\n\n")
        push!(sols, x.value)
    end

    close(f)
    redirect_stdout(originalSTDOUT)

    return sols
end 