# This file contains a library of predefined graphs. Based on the "Interesting" 4 node DAGs

# 2 parallel links
p2 = Graph([0 2;
            0 0])

# 3 parallel links
p3 = Graph([0 3;
            0 0])

# lollipop 
lol = Graph([0 1 0;
             0 0 2;
             0 0 0])

# augmented lollipop
a_lol = Graph([0 1 0;
               0 0 2;
               0 0 0])

# braess
braess = Graph([0 1 1 0;
                0 0 1 1;
                0 0 0 1;
                0 0 0 0])

# augmented braess
a_braess = Graph([0 1 1 1;
                0 0 1 1;
                0 0 0 1;
                0 0 0 0])

"""
Makes a square lattice graph of m * n with both incoming and outgoing 
edges between neighbours.
"""
function lattice(m,n)
    N = m*n # number of nodes
    A = zeros(Int, N,N) # start with zero matrix
    for k in 1:N
        (k%n != 0) && (A[k,k+1] = 1) # If not at right edge
        ((k-1)%n != 0) && (A[k, k-1] = 1) # If not at left edge
        (k+n <= N) && (A[k,k+n] = 1) # If not at bottom
        (k-n > 0) && (A[k,k-n] = 1) # If not at top
    end
    Graph(A)
end