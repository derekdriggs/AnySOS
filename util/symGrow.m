% Add redundant variables due to symmetry.
function [A,c,K] = symGrow(A,c,K)

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

projs = cell(1);

projs{1} = projMap(n);

P = blkdiag(projs{:});

A = [ A(:,1:nf+nl) A(:,nf+nl+1:end)*P];
            
c1 = c(1:nf+nl);
c2 = c(nf+nl+1:end);
c = [ c1; vec(P'*c2)];

end