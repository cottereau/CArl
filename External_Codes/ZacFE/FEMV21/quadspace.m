function x=quadspace(a,b,n)
%QUADSPACE  Quadratically spaced vector.
%   QUADSPACE(x1, x2) generates a row vector of 100 quadratically
%   linearly spaced points between x1 and x2.
%
%   QUADSPACE(x1, x2, N) generates N points between x1 and x2.
%   If N is negative the smallest space is at the enpoint x2,
%   otherwise it is at enpoint x1.
%
%   See also LINSPACE, LOGSPACE, LOGSP.

% Copyright (c) 2002-03-24, B. Rasmus Anthin.

error(nargchk(2,3,nargin))
if nargin<3,n=100;end
if n<-1
   n=-n;
   dx=cumsum(0:n-1);
   k=(b-a)/dx(n);
   x=fliplr(b-dx*k);
elseif n>1
   dx=cumsum(0:n-1);
   k=(b-a)/dx(n);
   x=a+dx*k;
else
   x=mean([a b]);
end