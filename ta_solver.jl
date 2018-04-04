using Convex, SCS, Gurobi
# Functions and calls to Convex.jl for solving Traffic Assignment

"""
    make_eq_constratints(rn::RoadNetwork, OD::AbstractMatrix, q::Float64, x::Variable)
    
Returns the equality constraints needed for the TA optimisation
of the RoadNetwork 'rn' and the demand level 'q'. The variable that is 
used in the problem must also be passed as an argument 'x' to the function

This version of the function uses sparse matrices and sparse arrays now...
"""
function make_eq_constratints(rn::RoadNetwork, OD::AbstractMatrix, q::Float64, x::Variable)
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
    make_ta_problem(rn::RoadNetwork, OD::AbstractMatrix, q::Float64, regime::String)
    
Generates convex optimisation problem from the graph g, the Origin-Destination matrix OD,
and the demand vector q. If only one OD pair, then q is a scalar. Regime is either "UE" (default)
which is "user equilibrium" or "SO": system optimal.

This should be a general method for multiple OD pairs, but it only works for 1 OD pair for now
This function uses Convex.jl
"""
function make_ta_problem(rn::RoadNetwork, OD, q::Float64, regime::String)
    
    a = rn.a
    b = rn.b
    m = num_edges(rn.g)

    #num_flows = sum(map(sign, rn.OD))
    x = Convex.Variable(m)
    
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
    ta_solve(rn::RoadNetwork, OD, q_range::Array{Float64,1}; regime="UE", solver=SCSSolver(verbose=false))
    
Returns solutions to the traffic assignment problem for a given range of demands
Calls function make_ta_problem

It returns a 2-dimensional array of the solution. The first row corresponds to the 
first element (first edge) and so on... There is a column for every demand step.

In theory it should work with any solver in MathProgBase, I have only tested SCS and Gurobi (which seems to 
be having some bugs though...)
"""
function ta_solve(rn::RoadNetwork, OD, q_range::Array{Float64,1}; regime="UE", solver=SCSSolver(verbose=false))
    
    println("Solving $regime STAP for d ∈ [$(q_range[1]), $(q_range[end])] ($(length(q_range)) step(s))\n")

    sols = Array{Float64,1}[]
    problem, x = make_ta_problem(rn, OD, q_range[1], regime)
    
    solve!(problem, solver)
    push!(sols, x.value[:])
    
    if length(q_range) == 1
        return sols[1][:]
    elseif length(q_range) > 1
        for q in q_range[2:end]
            problem.constraints[1] = make_eq_constratints(rn, OD, q, x)
            solve!(problem, warmstart=true)
            push!(sols, x.value[:])
        end
    end
    hcat(sols...)
end
ta_solve(rn::RoadNetwork, OD, q::Float64; regime="UE", solver=SCSSolver(verbose=false)) = ta_solve(rn, OD, [q], regime=regime, solver=solver)[:]


"""
    mixed_ta_solve(rn, od, demands::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-7, max_iters=50, warmstart_flag=true)

Solves a mixed equilibrium static traffic assignment. Where for the OD the demand is splilt into a proportion γ that attempts to
minimise total cost and (1 - γ) that tries to solve for user equilibrium. 

Defaults are: 
tolerance=1e-6
max_iters=50
"""
function mixed_ta_solve(rn, od, demands::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-7, max_iters=50, warmstart_flag=true)
    
    m = num_edges(rn.g)
    a = rn.a
    b = rn.b
    B = diagm(rn.b)    
    
    #Set up problems
    ue_d , so_d = demands
    x = Variable(m, Positive())
    y = Variable(m, Positive())
    
    ue_objective = a'*x + (2*B*y)'*x + 0.5*quadform(x, B)
    ue_problem = minimize(ue_objective, make_eq_constratints(rn, od, ue_d, x), x >= 0)

    so_objective = a'*y + quadform(y + x, B) + a'*x
    #so_objective = a'*y + quadform(y + x, B)
    so_problem = minimize(so_objective, make_eq_constratints(rn, od, so_d, y), y >= 0)
    
    x_init = ta_solve(rn, od, ue_d, regime="UE")
    y_init = ta_solve(rn, od, so_d, regime="SO")
    
    #x_init = zeros(m)
    #y_init = zeros(m)
    
    #First iteration
    fix!(x, x_init)
    solve!(so_problem, solver)
    free!(x)
    fix!(y, y_init)
    solve!(ue_problem, solver)
    free!(y)

    
    err = 10
    counter = 0    
    while (counter < max_iters) && (err > tolerance)
        
        old_x = copy(x.value[:])
        old_y = copy(y.value[:])
        
        fix!(y)
        solve!(ue_problem, solver, warmstart=warmstart_flag)
        free!(y)
        fix!(x)
        solve!(so_problem, solver, warmstart=warmstart_flag)
        free!(x)

        err = maximum([norm(x.value - old_x), norm(y.value - old_y)])
        counter += 1
    end
    #println("Iters:\t$(counter)\nError:\t$(err)\n")
    x.value[:], y.value[:], counter, err
end
"""
    mixed_ta_solve(rn, od, d::Float64, γ::Float64; solver=SCSSolver(verbose=false), tolerance=1e-6, max_iters=50, warmstart_flag=true)    
"""
function mixed_ta_solve(rn, od, d::Float64, γ::Float64; solver=SCSSolver(verbose=false), tolerance=1e-6, max_iters=50, warmstart_flag=true)
    ue_d = (1.0 - γ)*d
    so_d = γ*d             
    mixed_ta_solve(rn, od, [ue_d, so_d], solver=solver, tolerance=tolerance, max_iters=max_iters, warmstart_flag=warmstart_flag)
end

"""
    mixed_ta_solve(rn, od, d::Float64, γ_range::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-6, max_iters=50, warmstart_flag=true)    
"""
function mixed_ta_solve(rn, od, d::Float64, γ_range::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-6, max_iters=50, warmstart_flag=true)

    solution_flows_ue = Array{Float64,1}[]
    solution_flows_so = Array{Float64,1}[]
    solution_flows_agg = Array{Float64,1}[]
    iterations_to_converge = Int[]
    last_residual = Float64[]
     
    for γ in γ_range
        selfish_flows, altruistic_flows, iters, err = mixed_ta_solve(rn, od, d, γ, solver=solver, tolerance=tolerance, max_iters=max_iters)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
        push!(solution_flows_ue, selfish_flows)
        push!(solution_flows_so, altruistic_flows)
        push!(solution_flows_agg, selfish_flows + altruistic_flows)
        push!(iterations_to_converge, iters)
        push!(last_residual, err)

    end
    
    solution_flows_ue = unpack_sols(solution_flows_ue)
    solution_flows_so = unpack_sols(solution_flows_so) 
    solution_flows_agg = unpack_sols(solution_flows_agg)
    
    mixed_link_costs = [costs(rn, solution_flows_agg[:,i]) for i in 1:length(γ_range)]
    mixed_total_cost = [total_cost(rn, solution_flows_agg[:,i]) for i in 1:length(γ_range)]

    mixed_selfish_cost = [mixed_link_costs[i]⋅solution_flows_ue[:,i] for i in 1:length(γ_range)]
    mixed_alt_cost = [mixed_link_costs[i]⋅solution_flows_so[:,i] for i in 1:length(γ_range)]

    solution_flows_agg, solution_flows_so, solution_flows_ue, mixed_total_cost, mixed_selfish_cost, mixed_alt_cost, iterations_to_converge, last_residual
end
