% Remove redundant variables due to symmetry.
% This makes projection onto affine space faster.
% TODO: Assert dimensions of matrix.
function [A,c,K] = symShrink2(A,c,K)

n  = K.s;
nf = K.f;
nl = K.l;

if size(c,1) == 1 && size(c,2) >= 2
    error('c must be a column vector')
end

projMap = containers.Map('KeyType','uint32','ValueType','any');

for i = 1:1
    I = polymat2vec(reshape(1:n^2,n,n));
    projMap(n) = sparse((1:length(I))',I,ones(size(I)),length(I),n^2);
end

projs    = cell(1);
projs{1} = projMap(n);

P = blkdiag(projs{:});

ind = triu(ones(K.s)) == 1; % indices for matricising vectors

D = - 1/2*speye(K.s); % scale diagonal by 1/2
D = D(ind);
D = diag(D);
D = D + speye(size(D));

A = [ A(:,1:nf+nl) 2*A(:,nf+nl+1:end)*P'*D];
            
c1 = c(1:nf+nl);
c2 = c(nf+nl+1:end);
c = [ c1; vec(D'*P*c2)];

end