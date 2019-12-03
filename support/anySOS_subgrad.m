function [lambda,df] = anySOS_subgrad(x,e,K,useEig,varargin)
% Computes subgradient of the map
%
%   x -> max t   such that:   x - t*e is PSD
%
% using a generalised eigenvalue decomposition.
%
% Inputs:
%   x - evaluate subgradient at this point
%   e - strictly feasible point
%   K - structure containing number of variables in free, nonnegative,
%          and PSD cones in this order.
%   useEig - use eig (use lobpcg if false)
%   (optional) numEigs - block size in lobpcg
%   (optional) iter    - number of (cg) iterations in lobpcg

nf = K.f;
nn = K.l;
n  = K.s;

nvec = (n + 1) .* (n ./ 2);
nsum = 0;

if ~useEig
    numEigs = varargin{1};
    if numEigs < 0; numEigs = 5; end
    iter    = varargin{2};
    
    D    = zeros(numEigs*size(n,1),1);
    U    = zeros(numEigs*size(n,1));
    
    for i = 1:size(n,2)
        xs = full( my_smat( x(nf+nn+nsum+1:nf+nn+nsum+nvec(i)) ) );
        es = full( my_smat( e(nf+nn+nsum+1:nf+nn+nsum+nvec(i)) ) );

        [Ui,Di,fail] = lobpcg(randn(size(xs,1),numEigs),xs,es,1e-4,iter);

        D(numEigs*(i-1)+1:numEigs*i) = Di;
        U(sum(n(1:i-1))+1:sum(n(1:i)),numEigs*(i-1)+1:numEigs*i) = Ui;
        
        nsum = nsum + nvec(i);
    end

    D     = [x(nf+1:nf+nn); D];
    U     = [eye(n,nn), U];
 else
    
    D    = zeros(sum(n),1);
    U    = zeros(sum(n));
    
    for i = 1:size(n,2)
        xs = full( my_smat( x(nf+nn+nsum+1:nf+nn+nsum+nvec(i)) ) );
        es = full( my_smat( e(nf+nn+nsum+1:nf+nn+nsum+nvec(i)) ) );

        [Ui,Di]  = eig(xs,es);

        D(sum(n(1:i-1))+1:sum(n(1:i))) = diag(Di);
        U(sum(n(1:i-1))+1:sum(n(1:i)),sum(n(1:i-1))+1:sum(n(1:i))) = Ui;
        
        nsum = nsum + nvec(i);
    end

    D      = [x(nf+1:nf+nn); D];
    U      = [zeros(nn,size(U,2)+nn); zeros(size(U,1),nn), U];
    U(1:nn,1:nn) = eye(nn);
    
end

[lambda,index] = min(D);
evec      = U(:,index);

dfblk = evec*evec';
dfnn  = diag(dfblk(1:nn,1:nn));
dfblk(1:nn,:) = [];
dfblk(:,1:nn) = [];

df = [];

for i = 1:size(n,2)
    inds = sum(n(1:i-1))+1:sum(n(1:i));
    dfi = dfblk(inds,inds);
    df = [df; my_svec( dfi )];
end

df = [zeros(nf,1); dfnn; df];

end




