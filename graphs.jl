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

show(io::IO, n::Node) = print(io, "<$(n.index)> $(n.pos)")
Node() = Node(0, Float64[])


type Edge
    index::Int
    source::Node
    target::Node
    #params::Dict{AbstractString,Float64}
end

show(io::IO, e::Edge) = print(io, "<$(e.index)> ($(e.source.index) → $(e.target.index))")
#Edge(i::Int, s::Node, t::Node) = Edge(i::Int, s::Node, t::Node, Dict{AbstractString,Float64}())

length(e::Edge) = norm(e.target.pos - e.source.pos)



"""
Graph object (it has to be directed). Contains nodes, edges and dictionaries
of incoming and outgoing edges for each node.
"""
type Graph
    nodes::Array{Node,1}
    edges::Array{Edge,1}
    in_edges::Dict{Node, Array{Edge,1}} #This should probably change to make mem allocation more efficient....
    out_edges::Dict{Node, Array{Edge,1}}
end


"""
Show important information about graph, and unicode plot of the edges.
At the moment it assumes coordinates are known for the nodes... (this should not be assumed, and should be fixed soon...) 
"""
function show(io::IO, g::Graph)
    n = num_nodes(g)
    m = num_edges(g)
    output = "Graph:\nNodes - $(n) \nEdges - $(m)\n\n"
    print(io, output)
    
    if all(x -> length(x) == 2, Array{Float64,1}[g.nodes[i].pos for i in 1:n]) #If all nodes have 2D coordinates, then we can plot
       canvas = BrailleCanvas(60,25,
                              origin_x = 0.0, origin_y = 0.0,
                              width = 1.0, height = 1.0)

        xs, ys = node_positions(g)
        edge_coords = Tuple{Float64,Float64,Float64,Float64}[(g.edges[i].source.pos..., g.edges[i].target.pos...) for i in 1:m] 
        for i in 1:m
            lines!(canvas, edge_coords[i]...)
        end
        #points!(canvas, xs, ys, :red)
        print(io, canvas)
    end
end

"""
    Graph()
    
Makes an empty graph
"""
Graph() = Graph(Array{Node,1}[], Array{Edge,1}[], Dict{Node, Array{Edge,1}}(), Dict{Node, Array{Edge,1}}())


"""
    Graph(A::Array{Int64,2})
    
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

"""
    Graph(A::AbstractSparseMatrix)
    
Makes a graph from a sparse adjacency matrix.
"""
function Graph(A::AbstractSparseMatrix)
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
    
    nz_indexes = [ind2sub(A, i)  for i in find(A)]

    for idx in nz_indexes
        i = idx[1]
        j = idx[2]
        
        for k in 1:A[i,j]
            e = Edge(counter, nodes[i], nodes[j])
            push!(edges, e)
            push!(out_edges[nodes[i]], e)
            push!(in_edges[nodes[j]], e)
            counter += 1
        end
    end
    Graph(nodes, edges, in_edges, out_edges)
end



# Defining functions
#######################
#######################

"""
    add_node!(g::Graph, n::Node)
    
Adds node n to graph g.
"""
function add_node!(g::Graph, n::Node)
    indx = num_nodes(g) + 1
    n.index = indx

    push!(g.nodes, n)
    g.in_edges[n] = Edge[]
    g.out_edges[n] = Edge[] 
end

"""
    add_node!(g::Graph)
    
Adds a default node (with no position) to graph g (taking care of the indices).
"""
function add_node!(g::Graph)
    indx = num_nodes(g) + 1
    new_node =  Node(indx, Float64[])
    add_node!(g, new_node)
end

"""
    add_node!(g::Graph, pos::Array{Float64,1})
    
Adds a node to graph g at with position pos.
"""
function add_node!(g::Graph, pos::Array{Float64,1})
    indx = num_nodes(g) + 1
    n = Node(indx, pos)
    push!(g.nodes, n)
end

"""
    connect_net!(g::Graph, i::Int, j::Int)
    
Adds and edge pointing from node i to node j
"""
function connect_net!(g::Graph, i::Int, j::Int)
    m = num_edges(g)
    edge = Edge(m+1, g.nodes[i], g.nodes[j])
    push!(g.edges, edge)
    push!(g.out_edges[g.nodes[i]], edge)
    push!(g.in_edges[g.nodes[j]], edge)
end

"""
    add_edge!(g::Graph, e::Edge)
    
Adds an edge (given as Edge object) to the graph g. It also updates in_edges and out_edges.
"""
function add_edge!(g::Graph, e::Edge)
    m = num_edges(g)
    e.index = m + 1
    push!(g.edges, e)

    s = e.source
    t = e.target

    push!(g.out_edges[g.nodes[s.index]], e)
    push!(g.in_edges[g.nodes[t.index]], e)
end

"""
    add_edge!(g::Graph, i::Int, j::Int)
    
Adds a directed edge between given nodes. Where nodes are specified by their indices i and j.
"""
function add_edge!(g::Graph, i::Int, j::Int)
    connect_net!(g, i, j)
end

"""
    add_edge!(g::Graph, s::Node, t::Node)
    
Adds an edge between specified nodes. Nodes should belong to g already (probably...)
"""
function add_edge!(g::Graph, s::Node, t::Node)
    if s in g.nodes & t in g.nodes
        new_index = num_edges(g) + 1
        edge = Edge(new_index, s, t)
        push!(g.edges, e)
        push!(g.out_edges[s], e)
        push!(g.in_edges[t], e)
    else
        error("Nodes not in graph! (at least one of them...)")
    end
end

#"""
#Adds an edge to the graph that connects node i to node j,
#where i and j are the indices for 2 nodes in the graph g.
#
#This might actually not work if i remember correctly and doesnt it do
#the same as connect_net ?? (see if we need removing.)
#"""
#function connect!(g::Graph, i::Int, j::Int)
#    m = num_edges(g)
#    e = Edge(m+1, g.nodes[i], g.nodes[j])
#    add_edge!(g, e)
#end

"""
    num_nodes(g::Graph)
    
Returns the number of nodes of graph g.
"""
function num_nodes(g::Graph)
    length(g.nodes)
end

"""
    num_edges(g::Graph)
    
Returns the number of edges of graph g.
"""
function num_edges(g::Graph)
    length(g.edges)
end

"""
    in_edges_idx(n::Node, g::Graph)
    
Returns an array of the indices of the edges that are
incoming at node n.
"""
function in_edges_idx(n::Node, g::Graph)
    sort!([e.index for e in g.in_edges[n]])
end
in_edges_idx(i::Int, g::Graph) = in_edges_idx(g.nodes[i], g)

"""
    out_edges_idx(n::Node, g::Graph)
    
Returns an array of the indices of the edges that are
outgoing at node n.
"""
function out_edges_idx(n::Node, g::Graph)
    sort!([e.index for e in g.out_edges[n]])
end
out_edges_idx(i::Int, g::Graph) = out_edges_idx(g.nodes[i], g)

"""
    adjacency_matrix(g::Graph)
    
Returns the adjacency matrix of g. It is returned as a sparse matrix (SparseMatrixCSC).
For a full matrix, use adjacency_matrix_non_sparse.
"""
function adjacency_matrix(g::Graph)
    n = num_nodes(g)
    A = spzeros(n,n)
    for e in g.edges
        i, j = e.source.index, e.target.index
        A[i,j] += 1
    end
    A
end

"""
    adjacency_matrix_non_sparse(g::Graph)
    
Non-sparse version of the function adjacency matrix. 
Returns the adjacency matrix of the graph g. Returns a non-sparse matrix.
"""
function adjacency_matrix_non_sparse(g::Graph)
    A = adjacency_matrix(g)
    full(A)
end

"""
    incidence_matrix(g::Graph)
    
Returns the incidence matrix (SparseMatrixCSC) of 'g'. An n x m matrix where n is the number of nodes and m is the number of edges.
Convention for the incidence matrix:
M[i,j] = 1 if edge j is incoming at node i; 
M[i,j] = -1 if edge j is outgoing at node i
M[i,j] = 0 otherwiswe.
"""
function incidence_matrix(g::Graph)
    n = num_nodes(g)
    m = num_edges(g)
    M = spzeros(Int64, n, m)

    #Go through all nodes
    for i in 1:n
        out_e = out_edges_idx(i, g)
        in_e = in_edges_idx(i,g)
        #Outgoing edges
        for j in out_e
            M[i,j] = -1
        end
        #Incoming edges
        for j in in_e
            M[i,j] = 1
        end
    end
    M
end

"""
    incidence_matrix_non_sparse(g::Graph)
    
Non-sparse version of incidence_matrix.
"""
function incidence_matrix_non_sparse(g::Graph)
    full(incidence_matrix(g))
end

function node_positions(g::Graph)
    n = num_nodes(g)
    xs = Float64[g.nodes[i].pos[1] for i in 1:n]
    ys = Float64[g.nodes[i].pos[2] for i in 1:n]

    xs, ys
end
