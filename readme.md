# The Hierarchical Dirichlet Process

We sample the post of G0, G1, G2,...,GN explicitly.

The model is:

     G0  ~ DP(gamma, H)
     Gi  ~ DP(alpha, G0),      i = 1,...,N
     Mij ~ Multinomial(Gi),    j = 1,...,Ni
     Wij ~ Multinomial(T_Mij), j = 1,...,Ni
