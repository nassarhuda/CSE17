% [v,e] = DEGSEQ(p,dmax,dmin,n)
% Given a number of vertices, n, a maximum degree, dmax, minimum degree, dmin,
% and a power law exponent, p, this generates a degree sequence (contained in the length
% n vector v) satisfying a power law distribution with max degree d.

function [v,e] = degseq(p,dmax,dmin,n)
v = ones(n,1);
v = v.*dmin;

last = floor((dmin/dmax)^(-1/p));

for k = 1:last
	v(k) = dmax/(k^p);
end
