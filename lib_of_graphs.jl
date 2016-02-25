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
                0 0 0 0})

# augmented braess
braess = Graph([0 1 1 1;
                0 0 1 1;
                0 0 0 1;
                0 0 0 0])

