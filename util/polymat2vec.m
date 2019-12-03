function y=polymat2vec(x)
% function y=mss_s2v(x,z) from spotless toolbox.
%
% with nargin==1, re-arrange elements of matrix x into vector
%   x = [ x(1,1); x(1,2); x(2,2); x(1,3); x(2,3); x(3,3); ... x(m,m)],
% where m=min(size(x))
% otherwise


m=min(size(x));
x=x(1:m,1:m);
n=m*(m+1)/2;
ii=(1:n)';
b=ceil((sqrt(1+8*ii)-1)/2);
a=ii-b.*(b-1)/2;
y=x(a+m*(b-1));
if nargin>1,
    s=polymat2vec(polymat2vec(1:n));
    [s,jj]=sort(s);
    y=y(jj);
end
