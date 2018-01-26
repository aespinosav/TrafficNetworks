module TrafficNetworks

using Convex, SCS, Gadfly, DataFrames
import Base.show

export 
    Node, Edge, Graph, RoadNetwork,
    add_node!, add_edge!, connect!, #Do these *need* to be exported?
    num_nodes, num_edges, node_positions,
    in_edges_idx, out_edges_idx, #Do these *need* to be exported?
    adjacency_matrix, incidence_matrix,
    adjacency_matrix_non_sparse, incidence_matrix_non_sparse,
    replace_OD_matrix!, od_pairs, od_matrix_from_pair, make_ta_problem,
    ta_solve, ta_solve!,
    p2, p3, lol, a_lol, braess, a_braess, lattice, # These could be in an atlas type module by themselves
    total_cost, edge_costs, marginal_edge_costs,
    flows_data_frame, normflows_data_frame,
    plot_flows, plot_normalised_flows

include("graphs.jl")
include("road_networks.jl")
include("ta_solver.jl")
include("lib_of_graphs.jl")
include("analysis_functions.jl")
include("ta_plotting.jl")
end
