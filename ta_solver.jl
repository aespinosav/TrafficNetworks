# Functions and calls to Convex.jl for solving Traffic Assignment
"""
Returns the equality constraints needed for the TA optimisation
of the RoadNetwork 'rn' and the demand level 'q'. The variable that is 
used in the problem must also be passed as an argument 'x' to the function

This version of the function uses sparse matrices and sparse arrays now...
"""
function make_eq_constratints(rn::RoadNetwork, OD, q, x::Variable)
    n = num_nodes(rn)
    m = num_edges(rn)
    M = incidence_matrix(rn.g)

    d = sparsevec(zeros(n)) #changed this since it was giving problems when m != n... have to test but it looks right.
    flow_counter = 1
    
    indices = find(OD)
    for k in indices
        i, j = ind2sub(OD, k)
        if OD[i,j] > 0
            d[i] = -q[flow_counter]*OD[i,j]
            d[j] = q[flow_counter]*OD[i,j]
            flow_counter += 1
        end
    end
    eq_constraints = M*x == d
end

"""
Generates convex optimisation problem from the graph g, the Origin-Destination matrix OD,
and the demand vector q. If only one OD pair, then q is a scalar. Regime is either "UE" (default)
which is "user equilibrium" or "SO": system optimal.

This should be a general method for multiple OD pairs, but it only works for 1 OD pair for now

This function uses Convex.jl
"""
function make_ta_problem(rn::RoadNetwork, OD, q, regime)
    
    a = rn.a
    b = rn.b
    m = num_edges(rn.g)

    #num_flows = sum(map(sign, rn.OD))
    x = Variable(m)
    
    if regime == "UE"
        cost_function = dot(rn.a, x) + 0.5*quadform(x, diagm(rn.b))
    elseif regime == "SO"
        cost_function = dot(rn.a, x) + quadform(x, diagm(rn.b))
    else
        error("Regime must be either 'UE' or 'SO'...")
    end

    # Constraints to be passes to Convex.jl problem
    eq_constraints = make_eq_constratints(rn, OD, q, x)
    ineq_constraints = x >= 0

    #Make Convex problem
    problem = minimize(cost_function, eq_constraints, ineq_constraints)

    return problem, x
end

"""
Unpacks the sols output of ta_solve into a 2-dim array
"""
function unpack_sols(array_of_vectors)
    cat(2,array_of_vectors...)
end

"""
Returns solutions to the traffic assignment problem for a given range of demands: q_range.
Calls function make_ta_problem

It returns a 2-dimensional array of the solution. The first row corresponds to the 
first element (first edge) and so on... There is a column for every demand step.
"""
function ta_solve(rn::RoadNetwork, OD, q_range::Array{Float64,1}; regime="UE", logfile_name="log_ta_solve.txt")

    println("Will solve $regime, TA problem  for $(length(q_range)) values of demand...")

    sols = Array{Float64}[]

    #redirect output of Convex solver to a log file to avoid screen clutter
    #originalSTDOUT = STDOUT
    #f = open(logfile_name, "w")
    #redirect_stdout(f)
    println("Will solve $regime, TA problem  for $(length(q_range)) values of demand...\n")

    problem, x = make_ta_problem(rn, OD, q_range[1], regime)
    #first solution individually to start with the warmstart later
    solve!(problem)
    push!(sols, x.value)
    #Iterates next optimisation routines with warmstart
    if length(q_range) > 1
        for q in q_range[2:end]
            problem.constraints[1] = make_eq_constratints(rn, OD, q, x)
            solve!(problem, warmstart=true)
            push!(sols, x.value)
        end
    end

    # return stdout to original settings (closes logfile as well)
    #close(f)
    #redirect_stdout(originalSTDOUT)

    unpack_sols(sols)
end

ta_solve(rn::RoadNetwork, OD, q::Float64; regime="UE", logfile_name="log_ta_solve.txt") = ta_solve(rn, OD, [q]; regime=regime, logfile_name=logfile_name)[:]
