function A = create_graph(p,n,dmax,dmin)
% This function creates a graph of 
%       1. power law p, 
%       2. nb of nodes n,
%       3. max degree dmax,
%       4. min degree dmin

    trials = 5;
    rng(1);
    v = ceil(degseq(p,dmax,dmin,n));
    A = rand_graph_degree(v,trials);

end

