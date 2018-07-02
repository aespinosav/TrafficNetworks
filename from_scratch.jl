using TrafficNetworks, JuMP, Gurobi

#### Network definition ###

a = [0.1, 1.0, 0.05, 1.0, 0.1]
b = [1.0, 0.1, 0.05, 0.1, 1.0]

g = Graph([0 1 1 0;
           0 0 1 1;
           0 0 0 1;
           0 0 0 0])
           
r1 = [1.0,0,0,1,0]
r2 = [0.0,1,0,0,1]
r3 = [1.0,0,1,0,1]
routes = [r1, r2, r3]



#### Solve mixed STAP ####

m = num_edges(rn.g)
n = num_nodes(rn.g)
n_d = length(demand_range)
a = rn.a
b = rn.b
inc_mat = incidence_matrix(rn.g)

STAP_sel = Model(solver=GurobiSolver())
STAP_alt = Model(solver=GurobiSolver())

d_sel = 
s_alt
