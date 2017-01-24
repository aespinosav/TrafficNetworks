module TrafficNetworks

using Convex, SCS, Gadfly, DataFrames

import Base.show

export 
    Node, Edge, Graph, RoadNetwork,
    add_node!, add_edge!, connect!, num_nodes, num_edges,
    in_edges_idx, out_edges_idx, 
    adjacency_matrix, incidence_matrix,
    replace_OD_matrix!, make_ta_problem,
    ta_solve, ta_solve!, p2, p3, lol, a_lol, braess,
    a_braess, lattice, total_cost, edge_costs, marginal_edge_costs,
    flows_data_frame, normflows_data_frame,
    plot_flows, plot_normalised_flows, save_graph, load_graph

## Includes

include("graphs.jl")
include("road_networks.jl")
include("ta_solver.jl")
include("lib_of_graphs.jl")
include("analysis_functions.jl")
include("ta_plotting.jl")

end