#This part of the module utilises Convex which imports functions that interfere with JuMP so careful with that.
using Convex, SCS, Gurobi

"""
    make_eq_constratints(rn::RoadNetwork, OD::AbstractMatrix, q::Float64, x::Variable)
    
Returns the equality constraints needed for the TA optimisation
of the RoadNetwork 'rn' and the demand level 'q'. The variable that is 
used in the problem must also be passed as an argument 'x' to the function

This version of the function uses sparse matrices and sparse arrays now...

To extend for multiple origin-destination pairs the best way might be to create all the constraints
in one go. But the constraints interplay with each other (or do they) so it may require


           Important, here I might need to check for type instability to
           make sure it is not making things slow. With AbstractMatrix.
"""
function make_eq_constratints(rn::RoadNetwork, OD::AbstractMatrix, q::Float64, x::Variable)
    n = num_nodes(rn.g)
    m = num_edges(rn.g)
    M = incidence_matrix(rn.g)

    d = sparsevec(zeros(n)) #To change for multi od...
    flow_counter = 1 #For multi od
    
    indices = find(OD)
    for k in indices
        i, j = ind2sub(OD, k)
        if OD[i,j] > 0 #probably a faster way to check this?
            d[i] = -q[flow_counter]*OD[i,j]
            d[j] = q[flow_counter]*OD[i,j]
            flow_counter += 1 #For multi od
        end
    end
    eq_constraints = M*x == d
end

"""
    make_ta_problem(rn::RoadNetwork, OD::AbstractMatrix, q::Float64, regime::String)
    
Generates convex optimisation problem from the graph g, the Origin-Destination matrix OD,
and the demand vector q. If only one OD pair, then q is a scalar. Regime is either "UE" (default)
which is 'user equilibrium' or 'SO': system optimal.

This should be a general method for multiple OD pairs, but it only works for 1 OD pair for now
This function uses Convex.jl.

Maybe parametrizing the function by the type of matrix passed could make things much faster!
our street networks are sparse after all.
"""
function make_ta_problem(rn::RoadNetwork, OD, q::Float64, regime::String)
    
    a = rn.a
    b = rn.b
    m = num_edges(rn.g)
    x = Convex.Variable(m)

    # I should be able to use forward diff to calculate both UE and SO in one go
    # this would speed up since in my case we always want to compare both scenarios.
    # since they actually transform so does SO live in the tangent bundle???
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
    ta_solve(rn::RoadNetwork, OD, q_range::Array{Float64,1}; regime='UE', solver=SCSSolver(verbose=false))
    
Returns solutions to the traffic assignment problem for a given range of demands
Calls function make_ta_problem

It returns a 2-dimensional array of the solution. The first row corresponds to the 
first element (first edge) and so on... There is a column for every demand step.

In theory it should work with any solver in MathProgBase, I have only tested SCS and Gurobi (which seems to 
bee having some bugs though...)
"""
function ta_solve(rn::RoadNetwork, OD, q_range::Array{Float64,1}; regime="UE", solver=SCSSolver(verbose=false))
    
    m = length(rn.g.edges)
    n = length(rn.g.nodes)
    n_q = length(q_range) #number of demand steps (samples?)
    
    if n_q > 1
       println("Solving $regime STAP for d ∈ [$(q_range[1]), $(q_range[end])] ($(length(q_range)) step(s))\n")
    end
    
    sols = zeros(Float64, m, n_q)
    problem, x = make_ta_problem(rn, OD, q_range[1], regime)
    
    solve!(problem, solver)
    sols[:,1] = x.value[:]
    
    if length(q_range) > 1
        for (i,q) in enumerate(q_range[2:end])
            problem.constraints[1] = make_eq_constratints(rn, OD, q, x)
            solve!(problem, warmstart=true)
            sols[:,i+1] = x.value[:]
        end
    end
    sols
end
ta_solve(rn::RoadNetwork, OD, q::Float64; regime="UE", solver=SCSSolver(verbose=false)) = ta_solve(rn, OD, Float64[q], regime=regime, solver=solver)


"""
    mixed_ta_solve(rn, od, demands::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-7, max_iters=50, warmstart_flag=true)

Solves a mixed equilibrium static traffic assignment. Where for the OD the demand is splilt into a proportion γ that attempts to
minimise total cost and (1 - γ) that tries to solve for user equilibrium. 

Defaults are: 
tolerance=1e-6
max_iters=50
"""
function mixed_ta_solve(rn, od, demands::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-7, max_iters=100, warmstart_flag=true)
    
    m = num_edges(rn.g)
    a = rn.a
    b = rn.b
    B = diagm(rn.b)    
    
    #Set up problems
    ue_d , so_d = demands
    x = Variable(m, Positive())
    y = Variable(m, Positive())
    
    ue_objective = a'*x + (B*y)'*x + 0.5*quadform(x, B)
    #ue_objective = a'*x + sum(B*y.*x) + quadform(x, B)
    ue_problem = minimize(ue_objective, make_eq_constratints(rn, od, ue_d, x))#, x >= 0) # I think the Positive() argument takes care of this

    so_objective = a'*y + a'*x + quadform(y + x, B)
    #so_objective = a'*y + sum(2*b.*y.*x) + quadform(y, B)
    #so_objective = a'*y + quadform(y, B) + 2*(x')*B*y
    so_problem = minimize(so_objective, make_eq_constratints(rn, od, so_d, y))#, y >= 0) # I think the Positive() argument takes care of this
    
    x_init = ta_solve(rn, od, ue_d, regime="UE")
    y_init = ta_solve(rn, od, so_d, regime="SO")
    #x_init = zeros(m)
    #y_init = zeros(m)
    
    x.value = x_init
    y.value = y_init
    
    #First iteration
    fix!(x)
    solve!(so_problem, solver)
    free!(x)
    fix!(y)
    solve!(ue_problem, solver)
    free!(y)

    err = 10
    counter = 0    
    while (counter < max_iters) || (err > tolerance)
        
        old_x = x.value[:]
        old_y = y.value[:]
        
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
    mixed_ta_solve(rn, od, Float64[ue_d, so_d], solver=solver, tolerance=tolerance, max_iters=max_iters, warmstart_flag=warmstart_flag)
end

"""
    mixed_ta_solve(rn, od, d::Float64, γ_range::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-6, max_iters=50, warmstart_flag=true)
    
"""
function mixed_ta_solve(rn, od, d::Float64, γ_range::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-6, max_iters=50, warmstart_flag=true)
    m = num_edges(rn.g)
    a = rn.a
    b = rn.b
    n_γ = length(γ_range)

    sol_sel = zeros(Float64, m, n_γ)
    sol_alt = zeros(Float64, m, n_γ)
    its = zeros(Int, n_γ) #Iterations between so and ue flows
    last_change = zeros(Float64, n_γ) #Change in magnitude of sol vector in last iter

    for (i,γ) in enumerate(γ_range)
        sol_sel[:,i], sol_alt[:,i], its[i], last_change[i] =
        mixed_ta_solve(rn, od, d, γ, solver=solver, tolerance=tolerance, max_iters=max_iters, warmstart_flag=warmstart_flag)
    end
    sol_sel, sol_alt, its, last_change
end












function make_so_obj(y::Convex.Variable,  x, rn)
    a = rn.a
    B = diagm(rn.b)
    a'*y + quadform(y + x, B)
end
function make_ue_obj(x::Convex.Variable,  y, rn)
    a = rn.a
    B = diagm(rn.b)
    a'*x + (2*B*y)'*x + 0.5*quadform(x, B)
end

function alt_mixed_ta_solve(rn, od, demands::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-7, max_iters=100, warmstart_flag=true)
    
    m = num_edges(rn.g)
    a = rn.a
    b = rn.b
    B = diagm(rn.b)    
    
    #Set up problems
    ue_d , so_d = demands
    x = Variable(m, Positive())
    y = Variable(m, Positive())
    
    #ue_objective = a'*x + (2*B*y)'*x + 0.5*quadform(x, B)
    #ue_problem = minimize(ue_objective, make_eq_constratints(rn, od, ue_d, x))

    #so_objective = a'*y + quadform(y + x, B)
    #so_problem = minimize(so_objective, make_eq_constratints(rn, od, so_d, y))
    
    x_init = ta_solve(rn, od, ue_d, regime="UE")
    y_init = ta_solve(rn, od, so_d, regime="SO")
    #x_init = zeros(m)
    #y_init = zeros(m)    
    x.value = x_init
    y.value = y_init
    
    problem = minimize(make_so_obj(y, x.value[:], rn), make_eq_constratints(rn, od, so_d, y))
    solve!(problem, solver)
    problem.objective = make_ue_obj(x, y.value[:], rn)
    problem.constraints[1] = make_eq_constratints(rn, od, ue_d, x)
    solve!(problem, solver)

    err = 10
    counter = 0    
    while (counter < max_iters) || (err > tolerance)
        
        old_x = x.value[:]
        old_y = y.value[:]
        
        problem.objective = make_ue_obj(x, y.value[:], rn)
        problem.constraints[1] = make_eq_constratints(rn, od, ue_d, x)
        solve!(problem, solver)
        problem.objective = make_so_obj(y, x.value[:], rn)
        problem.constraints[1] = make_eq_constratints(rn, od, so_d, y)
        solve!(problem, solver)

        err = maximum([norm(x.value - old_x), norm(y.value - old_y)])
        counter += 1
    end
    x.value[:], y.value[:], counter, err
end


function alt_mixed_ta_solve(rn, od, d::Float64, γ::Float64; solver=SCSSolver(verbose=false), tolerance=1e-6, max_iters=100, warmstart_flag=true)
    ue_d = (1.0 - γ)*d
    so_d = γ*d             
    alt_mixed_ta_solve(rn, od, Float64[ue_d, so_d], solver=solver, tolerance=tolerance, max_iters=max_iters, warmstart_flag=warmstart_flag)
end


function alt_mixed_ta_solve(rn, od, d::Float64, γ_range::Array{Float64,1}; solver=SCSSolver(verbose=false), tolerance=1e-6, max_iters=100, warmstart_flag=true)
    m = num_edges(rn.g)
    a = rn.a
    b = rn.b
    n_γ = length(γ_range)

    sol_sel = zeros(Float64, m, n_γ)
    sol_alt = zeros(Float64, m, n_γ)
    its = zeros(Int, n_γ) #Iterations between so and ue flows
    last_change = zeros(Float64, n_γ) #Change in magnitude of sol vector in last iter

    for (i,γ) in enumerate(γ_range)
        sol_sel[:,i], sol_alt[:,i], its[i], last_change[i] =
        alt_mixed_ta_solve(rn, od, d, γ, solver=solver, tolerance=tolerance, max_iters=max_iters, warmstart_flag=warmstart_flag)
    end
    sol_sel, sol_alt, its, last_change
end
