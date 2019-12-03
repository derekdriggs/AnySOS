%%% TODO: check dimensions and assert
function Aflip = flip_A(A)
% To represent a linear operator acting on a symmetric matrix A_full(X),
%   it is common to remove entries redundant due to symmetry.
%   Some toolboxes (e.g. Spotless) represent X as a vector x1 of its
%   upper-triangular entries using column-major vectorisation.
%
%   The AnySOS convention is to scale the off-diagonal entries of X
%   by sqrt(2) and store it as a vector x2 of its lower-triangular entries 
%   (i.e. x2 = my_svec(X))
%
%   This function returns a matrix Aflip so that 
%   A(x1) = Aflip(x2) = A_full(X).


n     = (sqrt(1+size(A,2)*8)-1)/2;
ind   = triu(ones(n)) == 1;

temp = zeros(n,n);

x = 1:n*(n+1)/2;

% change indexing
temp(ind) = x;
temp = temp+triu(temp,1)';
reord = temp(ind')';

Aflip = A(:,reord);

% adjust scale
I = my_svec(speye(n));
Aflip(:,I ~= 1) = Aflip(:,I ~= 1)/sqrt(2);


end