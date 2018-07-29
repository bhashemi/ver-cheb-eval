function fx = ver_bary(x, fvals, xk, wk)
%VER_BARY   Extension of barycentric interpolation formula to interval
%           arithmetic.
% This is adapted from Chebfun.
%
% Input x is the evaluation point. 
% fvals are the values of the Chebyshev expansion at the barycentric nodes 
% xk, and wk are the barycentric weights.

% Parse inputs:
[n, m] = size(fvals);
sizex = size(x);
ndimsx = ndims(x);

if ( (m > 1) && (ndimsx > 2) )
    error('ver_bary:evalArrayAtNDArray', ...
        ['VER_BARY does not support evaluation of vectors of polynomials at ' ...
         'inputs with more than two dimensions.']);
end

% Default to Chebyshev nodes and barycentric weights:
if ( nargin < 3 )
    xk = verchebpts(n);
end

if ( nargin < 4 )
    [xk, wk] = verchebpts(n);
end

if ( ~all(sizex) )
    fx = x;
    return
end

% Check that input is a column vector:
if ( (ndimsx > 2) || (sizex(2) > 1) )
    x = x(:);
end

% The function is a constant.
if ( n == 1 )
    fx = repmat(fvals, length(x), 1);
    return
end

% The function is NaN.
if ( any(isnan(fvals)) )
    fx = NaN(length(x), m);
    return
end

% The main loop:
if ( length(x) < 4*length(xk) )  % Loop over evaluation points
    % Note: The value "4" here was detemined experimentally.

    % Initialise return value:
    fx = intval(zeros(size(x, 1), m));

    % Loop over x
    for j = 1:size(x,1)
        xx = wk ./ (x(j) - xk);
        fx(j,:) = (xx.'*fvals) / sum(xx);
    end
else                            % Loop over barycentric nodes
    % Initialise:
    num = intval(zeros(size(x, 1), m));
    denom = num;

    % Altenatively, loop over xk
    for j = 1:length(xk)
        tmp = (wk(j) ./ (x - xk(j)));
        num = num + tmp*fvals(j,:);
        denom = denom + tmp;
    end
    fx = num ./ denom;
end

% Try to clean up NaNs:
for k = find(isnan(fx(:,1)))'       % (Transpose as Matlab loops over columns)
    index = find(in(x(k),xk), 1);   % Warning: It is possible that x(k) 
                                    % does intersect with one of the 
                                    % barycentric nodes, but is NOT its 
                                    % subset. This happens e.g., if x(k) is
                                    % a widely perturbed chebpts itself. In
                                    % that case, I don't know any ways to
                                    % replace the resulting NaN entry with 
                                    % a verified correct value and we have 
                                    % to be happy with the NaN value in fx.
                                    % We can replace NaNs correcetly, only 
                                    % if x(k) is inside one of the nodes in
                                    % xk.
    if ( ~isempty(index) )
        fx(k,:) = fvals(index,:);   % Set entry/entries to the correct value
    end
end

% Try to clean up +-Infs caused by x = xk:
for k = find(isinf(fx(:,1)))'       % (Transpose as Matlab loops over columns)
    index = find(in(x(k),xk), 1);   % The same warning as above holds here 
                                    % too.
    if ( ~isempty(index) )
        fx(k,:) = fvals(index,:);   % Set entry/entries to the correct value
    end
end

% Reshape the output if possible:
if ( (m == 1) && ( (ndimsx > 2) || (sizex(2) > 1) ) )
    fx = reshape(fx, sizex);
elseif ( (m > 1) && ( (ndimsx == 2) || (sizex(2) > 1) ) )
    fx = reshape(fx, sizex(1), m*length(x)/sizex(1));
end

end
