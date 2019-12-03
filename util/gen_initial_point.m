% Create initial point satisfying affine constraints Ax = b
%  and c'*x = z.
function x0 = gen_initial_point(A,b,c,z)
    A(end+1,:) = c';
    b(end+1,:) = z;
   x0          = A\b;
end