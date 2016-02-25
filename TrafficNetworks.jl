module TrafficNetworks

using Convex, SCS

import Base.show

export 
    Node, Edge, Graph, RoadNetwork,
    add_node!, add_edge!, connect!, num_nodes, num_edges,
    in_edges_idx, out_edges_idx, 
    adjacency_matrix, incidence_matrix,
    replace_OD_matrix!, make_ta_problem,
    ta_solve, p2, p3, lol, a_lol, braess,
    a_braess

## Includes

include("graphs.jl")
include("road_networks.jl")
include("ta_solver.jl")

end