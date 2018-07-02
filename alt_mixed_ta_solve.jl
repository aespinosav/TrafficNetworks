using TrafficNetworks, Convex, SCS, Gurobi
import TrafficNetworks.mixed_ta_solve

function make_so_obj(y::Convex.Variable,  x)
    a'*y + quadform(y + x, B)
end
function make_ue_obj(x::Convex.Variable,  y)
    a'*x + (2*B*y)'*x + 0.5*quadform(x, B)
end

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
    #ue_problem = minimize(ue_objective, make_eq_constratints(rn, od, ue_d, x))

    so_objective = a'*y + quadform(y + x, B)
    #so_problem = minimize(so_objective, make_eq_constratints(rn, od, so_d, y))
    
    x_init = ta_solve(rn, od, ue_d, regime="UE")
    y_init = ta_solve(rn, od, so_d, regime="SO")
    #x_init = zeros(m)
    #y_init = zeros(m)    
    x.value = x_init
    y.value = y_init
    
    problem = minimize(make_so_obj(y, x.value[:]), make_eq_constratints(rn, od, so_d, y))
    solve!(problem, solver)
    problem.objective = make_ue_obj(x, y.value[:])
    problem.constraints[1] = make_eq_constratints(rn, od, ue_d, x)
    solve!(problem, solver)

    err = 10
    counter = 0    
    while (counter < max_iters) & (err > tolerance)
        
        old_x = x.value[:]
        old_y = y.value[:]
        
        problem.objective = make_ue_obj(x, y.value[:])
        problem.constraints[1] = make_eq_constratints(rn, od, ue_d, x)
        solve!(problem, solver)
        problem.objective = make_so_obj(y, x.value[:])
        problem.constraints[1] = make_eq_constratints(rn, od, so_d, y)
        solve!(problem, solver)

        err = maximum([norm(x.value - old_x), norm(y.value - old_y)])
        counter += 1
    end
    x.value[:], y.value[:], counter, err
end
