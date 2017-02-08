include("colnormout.jl")
include("splitRT.jl")
using Convex
using SCS
include("manage_procs.jl")

function aptrank_parallel(G::SparseMatrixCSC{Float64,Int64},
                Rtrain::SparseMatrixCSC{Int64,Int64},
                dH::SparseMatrixCSC{Float64,Int64},
                K::Int64,
                S::Int64,
                tao::Float64,
                diffusion_type::Int64,
                stat_choice::ASCIIString)

# AptRank parameters:
# K : total number of Markov Chain iterations
# S : number of shuffle time in splitting R into Rfit and Reval
# t : splitting percentage between Rfit vs. Reval
# diffusion_type:
# 1 : diffuse from G to H only; and
# 2 : diffuse between G and H.
# TODO: check types and check XXX comments

if ~isequal(stat_choice,"median") && ~isequal(stat_choice,"mean")
    error("Please choose a valid statistics option. Options are 'mean' or 'median'")
end

np = nprocs()
if np < 12
    addprocs(12-np)
end
np = nprocs()

m,n = size(Rtrain)
G = (G+G')/2
xc = zeros(Float64,K-1,S) # xc stores coefficients gamma in each shuffle
for s = 1:S
    Rfit,Reval = splitRT(Rtrain,tao)
    if diffusion_type == 1
        P = vcat(hcat(G,spzeros(m,n)),hcat(Rfit',dH))
    elseif diffusion_type == 2
        P = vcat(hcat(G,Rfit),hcat(Rfit',dH))
    else
        error("Diffusion type must be 1 or 2 way")
    end

    P = colnormout(P)
    Pt = P'
    @eval @everywhere Pt = $Pt
    P = 0
    gc()

    # Z = vcat(speye(m),spzeros(n,m))
    # algorithm requires computation of (1) Z = P*Z; and (2) Zh = Z[m+1:end,:]'
    # note that accessing columns is always faster than accessing rows
    # so I changed these two lines to the following:
    # Zt = Zt*Pt
    # and
    # Zh = Zt[:,m+1:end];

    # sample run:
    # julia> Z = vcat(speye(m),spzeros(n,m));
    # julia> Zt = hcat(speye(m),spzeros(m,n));
    # julia> P = colnormout(P);
    # julia> Pt = P';
    # julia> @time Z = P*Z;
    #    0.102851 seconds (25.59 k allocations: 48.152 MB, 11.02% gc time)
    # julia> @time Zt = Zt*Pt;
    #    0.149721 seconds (28.42 k allocations: 48.195 MB, 30.66% gc time)
    # julia> @time Zh = Z[m+1:end,:]';
    #    0.030649 seconds (21 allocations: 258.172 KB)
    # julia> @time Zh = Z[m+1:end,:]';
    #    0.028115 seconds (21 allocations: 258.172 KB)
    # julia> @time Zh = Zt[:,m+1:end];
    #    0.000083 seconds (13 allocations: 113.906 KB)
    # julia> @time Zh = Zt[:,m+1:end];
    #    0.000087 seconds (13 allocations: 113.906 KB)

    Zt = hcat(speye(m),spzeros(m,n))

    nb_rows = length(m+1:size(Pt,2))*size(Zt,1)
    A = zeros(Float64,nb_rows,K-1)

    N = size(Zt,1)
    bs = ceil(Int64,N/np)
    nblocks = ceil(Int64,N/bs)
    all_ranges = Array(UnitRange{Int64},nblocks)

    for i = 1:nblocks
        start = 1 + (i-1)*bs
        all_ranges[i] = start:min(start+bs-1,N)
    end

    for i = 1:np
        t = Zt[all_ranges[i],:]
        sendto(i,C=t)
    end


    for k = 1:K
        @everywhere C = C*Pt

        k == 1 && continue

        Zh = spzeros(size(Zt,1),length(m+1:size(Zt,2)))
        for i = 1:np
            Ci = getfrom(i, :C)
            Zh[all_ranges[i],:] = Ci[:,m+1:size(Zt,2)]
        end

        ii,jj,vv = findnz(Zh)
        Zh_size = size(Zh,1)
        rowids = ii+(jj-1)*Zh_size

        # the three lines of code above is due to Zh[:] being too slow
        # another way to do this is,
        # vv = reshape(Zh,prod(size(Zh)),1)
        # A[:,k-1] = vv
        # but with multiple time tests the two methods performed similarly
        # with slight preference for the method I ended up choosing
        # the preference was due to A being already constructed and
        # is sparse, so we really just need a small number of the elements in Zh

        A[rowids,k-1] = vv
        Zh = 0
        gc()
    end

    Qa,Ra = qr(A)
    A = 0
    gc()
    d = size(Ra,2)

    x = Variable(d)

    # the line of code below is a replacement to b = Reval[:] since [:] is too slow
    b = reshape(Reval,prod(size(Reval)),1)

    problem = minimize(norm2(Ra*x - Qa'*b), x >= 0, sum(x)==1)
    solve!(problem)
    xc[:,s] = vec(x.value)
end

# We've found values for x, now we need to solve the PageRank problem

    if isequal(stat_choice,"median")
        xopt = median(xc,2)
    elseif isequal(stat_choice,"mean")
        xopt = mean(xc,2)
    end

    xopt = xopt/sum(xopt)

    if diffusion_type == 1
        P = vcat(hcat(G,spzeros(m,n)),hcat(Rtrain',dH))
    elseif diffusion_type == 2
        P = vcat(hcat(G,Rtrain),hcat(Rtrain',dH))
    else
        error("Diffusion type must be 1 or 2 way")
    end

    P = colnormout(P)
    Pt = P'
    @eval @everywhere Pt = $Pt
    P = 0
    gc()

    Zt = hcat(speye(m),spzeros(m,n))

    nb_rows = length(m+1:size(Pt,2))*size(Zt,1)
    A = zeros(Float64,nb_rows,K-1)

    N = size(Zt,1)
    bs = ceil(Int64,N/np)
    nblocks = ceil(Int64,N/bs)
    all_ranges = Array(UnitRange{Int64},nblocks)

    for i = 1:nblocks
        start = 1 + (i-1)*bs
        all_ranges[i] = start:min(start+bs-1,N)
    end

    for i = 1:np
        t = Zt[all_ranges[i],:]
        sendto(i,C=t)
    end

    for k = 1:K
        @everywhere C = C*Pt

        k == 1 && continue

        Zh = spzeros(size(Zt,1),length(m+1:size(Zt,2)))
        for i = 1:np
            Ci = getfrom(i, :C)
            Zh[all_ranges[i],:] = Ci[:,m+1:size(Zt,2)]
        end

        ii,jj,vv = findnz(Zh)
        Zh_size = size(Zh,1)
        rowids = ii+(jj-1)*Zh_size

        A[rowids,k-1] = vv
        Zh = 0
        gc()
    end

    Xa = A*xopt
    Xa = reshape(Xa,m,n)

    A = 0
    gc()

    return (xc,Xa)
end
