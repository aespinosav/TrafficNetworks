#import Base.show

# Defining Types
######################
######################

type Node
    index::Int
    pos::Array{Float64,1}
    #in_edges::Edge[]
    #out_edges::Edge[]
end
# Default empty node
Node() = Node(0, Float64[])

type Edge
    index::Int
    source::Node
    target::Node
    #params::Dict{AbstractString,Float64}
end
#Edge(i::Int, s::Node, t::Node) = Edge(i::Int, s::Node, t::Node, Dict{AbstractString,Float64}())

"""
Graph object (it has to be directed). Contains nodes, edges and dictionaries
of incoming and outgoing edges for each node.
"""
type Graph
    nodes::Array{Node,1}
    edges::Array{Edge,1}
    in_edges::Dict{Node, Array{Edge,1}}
    out_edges::Dict{Node, Array{Edge,1}}
end

"""
Makes an empty graph
"""
Graph() = Graph(Array{Node,1}[], Array{Edge,1}[], Dict{Node, Array{Edge,1}}(), Dict{Node, Array{Edge,1}}())

"""
Makes a (multi di) graph from an adjacency matrix. The weights in the adjacency matrix should be 
integers and correspond to the number of edges between the given pair of nodes. 
"""
function Graph(A::Array{Int64,2})
    if size(A)[1] != size(A)[2]
        error("Adjacency matrix is not square")
    end

    in_edges = Dict{Node, Array{Edge,1}}()
    out_edges = Dict{Node, Array{Edge,1}}()

    number_of_nodes = size(A)[1]
    nodes = Node[]
    for i in 1:number_of_nodes
        push!(nodes, Node(i,Float64[]))
        out_edges[nodes[i]] = Edge[]
        in_edges[nodes[i]] = Edge[]
    end

    edges = Edge[]
    counter = 1
    for i in 1:number_of_nodes
        for j in 1:number_of_nodes
            for k in 1:A[i,j]

                e = Edge(counter, nodes[i], nodes[j])
                push!(edges, e)
                push!(out_edges[nodes[i]], e)
                push!(in_edges[nodes[j]], e)
                counter += 1
            end
        end
    end
    Graph(nodes, edges, in_edges, out_edges)
end






# show methods for neatness in REPL
show(io::IO, n::Node) = print(io, "<$(n.index)> $(n.pos)")
show(io::IO, e::Edge) = print(io, "<$(e.index)>($(e.source.index) → $(e.target.index))")
function show(io::IO, g::Graph)
    output = "Graph:\nNodes - $(num_nodes(g)) \nEdges - $(num_edges(g))"
    print(io, output)
end

# Defining functions
#######################
#######################

"""
Adds node n to graph g.
"""
function add_node!(g::Graph, n::Node)
    indx = num_nodes(g) + 1
    n.index = indx

    push!(g.nodes, n)
    g.in_edges[n] = Edge[]
    g.out_edges[n] = Edge[]
    return
end

"""
Adds a default node (with no position) to graph g (taking care of the indices).
"""
function add_node!(g::Graph)
    indx = num_nodes(g) + 1
    push!(g.nodes, Node(indx, Float64[]))
end

"""
Adds a node to graph g at with position pos.
"""
function add_node!(g::Graph, pos::Array{Float64,1})
    indx = num_nodes(g) + 1
    n = Node(indx, pos)
    push!(g.nodes, n)
end

"""
Adds an edge to the graph g. It actually reconstructs the graph
to get the edge indexes to follow the standard ordering.

This way of doing it might not be a good idea...
"""
function add_edge!(g::Graph, e::Edge)
    m = num_edges(g)
    e.index = m + 1
    push!(g.edges, e)
    A = adjacency_matrix(g)
    g = Graph(A)
end

"""
Adds an edge to the graph that connects node i to node j,
where i and j are the indices for 2 nodes in the graph g.
"""
function connect!(g::Graph, i::Int, j::Int)
    m = num_edges(g)
    e = Edge(m+1, g.nodes[i], g.nodes[j])
    add_edge!(g, e)
end



"""
Returns the number of nodes of graph g.
"""
function num_nodes(g::Graph)
    length(g.nodes)
end

"""
Returns the number of edges of graph g.
"""
function num_edges(g::Graph)
    length(g.edges)
end

"""
Returns an array of the indices of the edges that are
incoming at node n.
"""
function in_edges_idx(n::Node, g::Graph)
    sort!([e.index for e in g.in_edges[n]])
end
in_edges_idx(i::Int, g::Graph) = in_edges_idx(g.nodes[i], g)

"""
Returns an array of the indices of the edges that are
outgoing at node n.
"""
function out_edges_idx(n::Node, g::Graph)
    sort!([e.index for e in g.out_edges[n]])
end
out_edges_idx(i::Int, g::Graph) = out_edges_idx(g.nodes[i], g)


"""
Returns the adjacency matrix of the graph g.
"""
function adjacency_matrix(g::Graph)
    A = zeros(Int, num_nodes(g), num_nodes(g))
    for e in g.edges
        A[e.source.index, e.target.index] += 1
    end
    A
end


"""
Returns the incidence matrix of g. An n x m matrix where n is the nummber of 
nodes and m is the number of edges.

Convention for the incidence matrix:
M[i,j] = 1 if edge j is incoming at node i; 
M[i,j] = -1 if edge j is outgoing at node i
M[i,j] = 0 otherwiswe.
"""
function incidence_matrix(g::Graph)

    M = zeros(Int, num_nodes(g),num_edges(g))

    for i in 1:num_nodes(g)
        out_e = out_edges_idx(i, g)
        for j in out_e
            M[i,j] = -1
        end
    end

    for i in 1:num_nodes(g)
        in_e = in_edges_idx(i, g)
        for j in in_e
            M[i,j] = 1
        end
    end
    
    return M
end

"""
Saves a graph to a textfile that can be loaded with load_graph.
To work properly the nodes should have positions.
Tildes "~" indicate beggining and end of node and edge sections
"""
function save_graph(g, filename)

    header = "#TrafficNetworks graph (can be loaded into julia with load_graph)\n"
    node_header = "#Nodes\n#Index\tx\ty\n"
    edge_header = "#Edges\n#source\ttarget\n"

    indices_nodes = [n.index for n in g.nodes]
    xs = [n.pos[1] for n in g.nodes]
    ys = [n.pos[2] for n in g.nodes]

    indices_edges = [e.index for e in g.edges]
    sources = [e.source.index for e in g.edges]
    targets = [e.target.index for e in g.edges]

    open(filename, "w") do f
        write(f, header)
        write(f, node_header)
        write(f, "~\n")
        for i in 1:length(g.nodes)
            write(f, "$(indices_nodes[i])\t$(xs[i])\t$(ys[i])\n")
        end
        write(f, "~\n")
        write(f, edge_header)
        write(f, "~\n")
        for i in 1:length(g.edges)
            write(f, "$(indices_edges[i])\t$(sources[i])\t$(targets[i])\n")
        end
        write(f, "~\n")
    end
end

"""
Loads a graph file and returns a
"""
function load_graph(filename)

    g = Graph()

    ln = "go"
    open(filename, "r") do f
        while ln[1] != '~'
            ln = readline(f)
        end
        #Node adding start
        ln = readline(f)
        while ln[1] != '~'
            line = split(ln)
            ind = parse(Int, line[1])
            x = float(line[2])
            y = float(line[3])

            add_node!(g, Node(ind, [x, y]))

            ln = readline(f)
        end
        ln = readline(f)
        ln = readline(f)
        ln = readline(f)
        ln = readline(f)
        while ln[1] != '~'
            line = split(ln)
            i = parse(Int, line[2])
            j = parse(Int, line[3])

            connect!(g, i, j)
            ln = readline(f)
        end
    end
    return g
end