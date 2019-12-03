function t = anySOS_line_search(x,e,useEig,K)
% Solves the problem
%
%   max t   such that:   x - t*e is PSD
%
% using a generalised eigenvalue decomposition.

nf = K.f;
nn = K.l;
n  = K.s;

nvec = (n + 1) .* (n ./ 2);
nsum = 0;

ti = zeros(size(n));

if useEig
    for i = 1:size(n,2)
        xs = my_smat( x(nf+nn+nsum+1:nf+nn+nsum+nvec(i)) );
        es = my_smat( e(nf+nn+nsum+1:nf+nn+nsum+nvec(i)) );

        ti(i)  = min(eig(xs,es));
        nsum   = nsum + nvec(i);
    end
else
    for i = 1:size(n,2)
    
        nsum          = nsum + n(i);

        xs      = my_smat( x(nf+nn+nsum+1:nsum+n(i)) );
        es      = my_smat( e(nf+nn+nsum+1:nsum+n(i)) );
        
        [~,D,fail] = lobpcg(randn(size(x,1),5),x,e,1e-4,50);
        ti = min(D);
    end
end

t = min(ti);

end
