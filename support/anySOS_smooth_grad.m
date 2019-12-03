function [lambda,df] = anySOS_smooth_grad(x,e,mu,K,varargin)
% Computes subgradient of the map
%
%   x -> mu log( exp( \sum_i lambda_i(x)/mu ) ),
%
% using a generalised eigenvalue decomposition. We use a first-order 
% approximation for numerically stability.
%
% Inputs:
%   x  - evaluate subgradient at this point
%   e  - strictly feasible point
%   mu - smoothing parameter
%   K  - structure containing number of variables in free, nonnegative,
%          and PSD cones in this order.
%   useEig - use eig (use lobpcg if false)
%   (optional) numEigs - block size in lobpcg
%   (optional) iter    - number of (cg) iterations in lobpcg

nf = K.f;
nn = K.l;
n  = K.s;

nvec = (n + 1) .* (n ./ 2);
nsum = 0;

D    = zeros(sum(n),1);
U    = zeros(sum(n));


if nargin > 4
    numEigs = varargin{1};
    iter    = varargin{2};
    [U,D,fail] = lobpcg(randn(size(x,1),numEigs),x,e,1e-4,iter);
    D     = [x(nf+1:nf+nn); D];
    U     = [eye(n,nn), U];
    lambda     = min(D);
    D          = D - max(D); % Shift for numerical stability.
    D          = -D./mu;
else
    for i = 1:size(n,2)
        xs = my_smat( x(nf+nn+nsum+1:nf+nn+nsum+nvec(i)) );
        es = my_smat( e(nf+nn+nsum+1:nf+nn+nsum+nvec(i)) );
        
        [Ui,Di]  = eig(xs,es);

        D(sum(n(1:i-1))+1:sum(n(1:i))) = diag(Di);
        U(sum(n(1:i-1))+1:sum(n(1:i)),sum(n(1:i-1))+1:sum(n(1:i))) = Ui;
        
        nsum = nsum + nvec(i);
    end

    D      = [x(nf+1:nf+nn); D];
    U      = [zeros(nn,size(U,2)+nn); zeros(size(U,1),nn), U];
    U(1:nn,1:nn) = eye(nn);
    lambda = min(D);
    D      = -D./mu;
    
end

% TODO: Keep U in blocks

[D,index] = sort(D,'descend');
D         = D - D(1);           % For numerical stability
U         = U(:,index);
ExpD      = exp(D);
D_new     = ExpD./sum(ExpD);

dfblk = U*diag(D_new)*U';
dfnn  = diag(dfblk(1:nn,1:nn));
dfblk(1:nn,:) = [];
dfblk(:,1:nn) = [];

df = [];

for i = 1:size(n,2)
    index = sum(n(1:i-1))+1:sum(n(1:i));
    dfi = dfblk(index,index);
    df = [df; my_svec( dfi )];
end
df = [zeros(nf,1); dfnn; df];
end







