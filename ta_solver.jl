# Functions and calls to Convex.jl for solving Traffic Assignment

"""
Returns the equality constraints needed for the TA optimisation
of the RoadNetwork rn and the demand level q. The variable that is 
used in the problem must also be passed as an argument x to the function
"""
function make_eq_constratints(rn::RoadNetwork, q, x::Variable)
    n = num_nodes(rn)
    m = num_edges(rn)
    M = incidence_matrix(rn.g)

    d = zeros(n) #changed this since it was giving problems when m != n... have to tes but it looks right.
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

    eq_constraints = incidence_matrix(rn.g)*x == d
end

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
    eq_constraints = make_eq_constratints(rn, q, x)
    ineq_constraints = x >= 0

    problem = minimize(cost_function, eq_constraints, ineq_constraints)

    return problem, x
end

"""
Returns solutions to the traffic assignment problem for a given range of demands: q_range.
Calls function make_ta_problem

It returns a 2-dimensional array of the solution. The first row corresponds to the 
first element (first edge) and so on... There is a column for every demand step.
"""
function ta_solve!(ta::TrafficAssignment; regime="UE", logfile_name="log_ta_solve.txt")
    
    println("Will solve $regime, TA problem  for $(length(ta.demand_range)) values of demand...\n")
    
    rn = ta.rn
    ta.regime = regime
    
    #redirect output of Convex solver to a log file to avoid screen clutter
    originalSTDOUT = STDOUT
    f = open(logfile_name, "w")
    redirect_stdout(f)
    println("Will solve $regime, TA problem  for $(length(ta.demand_range)) values of demand...\n")
    
    # make optimisation problem (Convex.jl)
    problem, x = make_ta_problem(rn, ta.demand_range[1], regime)
    
    #first solution
    solve!(problem)
    row = [ta.demand_range[1], x.value...]
    push!(ta.flows, row)
    
    # iterates next optimisation routines with warmstart
    if length(ta.demand_range) > 1
        for q in ta.demand_range[2:end]
            problem.constraints[1] = make_eq_constratints(rn, q, x)
            solve!(problem, warmstart=true)
            row = [q, x.value...]
            push!(ta.flows, row)
        end
    end

    # return stdout to original settings (closes logfile as well)
    close(f)
    redirect_stdout(originalSTDOUT)
end

#ta_solve(rn::RoadNetwork, q::Float64; regime="UE", logfile_name="log_ta_solve.txt") = ta_solve(rn, [q], regime=regime, logfile_name=logfile_name)


# """
# Returns solutions to the traffic assignment problem for a given range of demands: q_range.
# Calls function make_ta_problem
# 
# It returns a 2-dimensional array of the solution. The first row corresponds to the 
# first element (first edge) and so on... There is a column for every demand step.
# 
# This function also updates the RoadNetwork rn to contain the demand range and 
# the flow solutions.
# 
# ta_solve(rn, q_range, regime, logfile_name)
#                      |  
#                       optional from here ->
# """
# function ta_solve!(rn::RoadNetwork, q_range::Array{Float64,1}; regime="UE", logfile_name="log_ta_solve.txt")
#     println("Will solve $regime, TA problem  for $(length(q_range)) values of demand...\n")
#     
#     rn.demand_range = q_range
# 
#     sols = Array{Float64}[]
# 
#     #redirect output of Convex solver to a log file to avoid screen clutter
#     originalSTDOUT = STDOUT
#     f = open(logfile_name, "w")
#     redirect_stdout(f)
# 
#     problem, x = make_ta_problem(rn, q_range[1], regime)
#     #first solution
#     solve!(problem)
#     push!(sols, x.value)
#     # iterates next optimisation routines with warmstart
#     if length(q_range) > 1
#         for q in q_range[2:end] # We are omitting the first value in q_Range as it is assumed you are starting from 0. This should probably change...
#             problem.constraints[1] = make_eq_constratints(rn, q, x)
#             solve!(problem, warmstart=true)
#             push!(sols, x.value)
#         end
#     end
# 
#     # return stdout to original settings (closes logfile as well)
#     close(f)
#     redirect_stdout(originalSTDOUT)
# 
#     flows = unpack_sols(sols)
#     if regime=="UE"
#         rn.flows_ue = flows
#     else
#         rn.flows_so = flows
#     end
# 
#     return flows
# end 
# 
# ta_solve!(rn::RoadNetwork, q::Float64; regime="UE", logfile_name="log_ta_solve.txt") = ta_solve!(rn, [q], regime=regime, logfile_name=logfile_name)
# 
# """
# Unpacks the sols output of ta_solve into a 2-dim array
# """
# function unpack_sols(array_of_vectors)
#     cat(2,array_of_vectors...)
# end