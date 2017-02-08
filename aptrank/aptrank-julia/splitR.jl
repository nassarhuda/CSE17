function splitR(R,rho)

# R is a sparse matrix and t is a value in (0,1)
# # TODO: error checking R must be sparse, t must be between (0,1)
# # TODO: add types
# % Split R into Rtrain and Rtest; rho, spliting percentage in [0,1]
# if rho < 0 || rho > 1
#     error('rho must be in [0,1]!');
# end
#
# if ~issparse(R)
#     R = sparse(R);
# end

m,n = size(R)
ei,ej,ev = findnz(R)
len = length(ev)
seed = time()
r = MersenneTwister(round(Int64,seed))
a = randperm(r,len)
nz = floor(Int,rho*len)
p = a[1:nz]
cp = setdiff(collect(1:len),p);

Rtrain = sparse(ei[p],ej[p],ev[p],m,n)
Rtest = sparse(ei[cp],ej[cp],ev[cp],m,n)

#assert(isequal(Rtrain+Rtest,R))

return(Rtrain,Rtest)
end
