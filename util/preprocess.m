% TODO: Make x0 a cell if provided
function [A,b,c,K,e,x0] = preprocess(A,b,c,K,e,options)

% Make sure we have column vectors
if size(b,1) == 1 && size(b,2) > 1
    b = b';
end

if size(c,1) == 1 && size(c,2) > 1
    c = c';
end


% Sparsify if the data are sparse
[n,m] = size(A);
densityTol = 0.4;
if nnz(A)/(m*n) > densityTol
    A = full(A);
else
    A = sparse(A);
end
if nnz(b)/m > densityTol
    b = full(b(:));
else
    b = sparse(b(:));
end
if nnz(c)/n > densityTol
    c = full(c);
else
    c = sparse(c);
end


% check that only free, zero, non-negative, quadratic cone and SDP variables
% are included
if isempty(K); error('Empty cone constraint types.'); end
if(~all(ismember(fieldnames(K),{'f','l','q','s','r','scomplex','ycomplex'})))
    error('Unsupported cone constraint types.');
end

if isfield(K,'r') || isfield(K,'q')
    if isempty(K.r); K.r = 0; end
    if isempty(K.q); K.q = 0; end
    if K.r ~= 0 || K.q ~= 0
        error('Cannot handle Lonentz cone constraints at this time.')
    end
end

if isfield(K,'scomplex') || isfield(K,'ycomplex')
    if isempty(K.scomplex); K.scomplex = 0; end
    if isempty(K.rcomplex); K.rcomplex = 0; end
    if K.scomplex ~= 0 || K.rcomplex ~= 0
        error('Cannot complex cone constraints at this time.')
    end
end

if ~isfield(K,'f') || isempty(K.f); K.f = 0; end
if ~isfield(K,'l') || isempty(K.l); K.l = 0; end
if ~isfield(K,'s') || isempty(K.s); K.s = 0; end


% Check if cone dimensions match with data
assert( size(A,1) == size(b,1), 'Dimension mismatch in SDP data');
assert( size(A,2) == size(c,1), 'Dimension mismatch in SDP data');

% Be sure that e satisfies affine constraint
assert( norm( A*e - b ) <= 1e-6, 'Matrix e does not satisfy affine constraint');

% Hamza: Check that e is strictly feasible?
assert( all(e((K.f+1):(K.f+K.l)) > 0) , 'e is not strictly feasible');
top = K.f+K.l;
for i=1:length(K.s)
    si = K.s(i);
    ind = (top+1):(top+si*(si+1)/2);
    assert( min(eig(my_smat( e(ind) ))) > 0 , 'e is not strictly feasible');
    top = top+si*(si+1)/2;
end
assert(top == length(e));

if ~isfield(options,'z') || isempty(options.z)
    z = c'*e - abs(c'*e);
else
    z = options.z;
end

if ~isfield(options,'x0') % Check if we need to generate initial point
    x0 = gen_initial_point(A,b,c,z);
else
    x0 = options.x0;
end

end
