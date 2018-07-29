function [x, baryW] = verchebpts(n, dom)
% VERCHEBPTS(N) returns N intervals containing Chebyshev points of the 2nd 
% kind in [-1,1].
% Adapated from Chebfun.

if ( n == 0 )     % Special case (no points)
    x = []; 
    
elseif ( n == 1 ) % Special case (single point)
    x = intval(0); 
    
else              % General case
    % Chebyshev points:
    m = n - 1;
    x = sin((intval('pi')*(-m:2:m))/(2*m)).';  
    % (Use of sine enforces symmetry.)    
end

if nargin > 1
    % Scale the nodes:
    x = dom(2)*(x + 1)/2 + dom(1)*(1 - x)/2;
end

if nargout > 1
% The following comes from chebtech2.barywts(n).

    if ( n == 0 )                      % Special case (no points)
        baryW = [];
    elseif ( n == 1 )                  % Special case (single point)
        baryW = 1;
    else
        % General case
        baryW = [ones(n-1,1); 1/2];    % Note v(end) is positive.
        baryW(end-1:-2:1) = -1; 
        baryW(1) = baryW(1)/2;
        
    end
end

end