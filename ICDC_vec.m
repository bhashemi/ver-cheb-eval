function [z, bInt] = ICDC_vec(t, coeffs, prec)

%% Step 1: Floating point approximation to ptilde = p(t) using Clenshaw: 
% approximately solve A x b = coeffs
[ptilde, b] = clenshaw(mid(t), mid(coeffs)); % Compute b s.t. A*b = coeffs,
% without ptilde and b should both be POINT data, not intervals!

%% Step 2: Enclose the residual r = coeffs - A*b without explicit A
if strcmp(prec, 'double')
    r = resEnc(t, b, coeffs, prec);
elseif strcmp(prec, 'quad')
    for i = 1:size(t,1)
        r(:,i) = resEnc(t(i), b(:,i), coeffs, prec);
    end
end

%% Step 3: compute the interval solution with the residual as RHS    
[y, bInt] = clenshaw(intval(t), r);

%% Step 4: Add the enclosure for the defect to the approximate solution
z = ptilde + y;

end

function [p, b] = clenshaw(x, coeffs)
% Clenshaw's algorithm to evaluate a polynomial p(x) at x in [-1, +1]
% If x is a scalar, then p = p(x) is a scalar as well and the recurrence 
% history b is a vector.
% x can be a vector in which case b is a matrix.

bk1 = 0*x;  % same type as x
bk2 = bk1;  % same type as x
x2 = 2*x;
len = size(coeffs,1); % degree of poly
b = zeros(len,size(x,1));
if isintval(x)
    b = intval(b);
end

for k = len:-1:1
    bk = coeffs(k) + (x2.*bk1 - bk2);
    bk2 = bk1; 
    bk1 = bk;
    b(k,:) = bk;
end
p = b(1,:) - b(2,:).*x';
end

%%
function r = resEnc(t, b, coeffs, prec)
n = size(coeffs,1);
m = size(t,1);

if strcmp(prec, 'double')
    b = intval(b);
    r = intval(zeros(n,m));
    for i=n-2:-1:1
        r(i,:) = ((coeffs(i) - b(i,:)) + 2*(intval(t).'.*b(i+1,:))) - b(i+2,:);
        % The order of summation above does matter in the quality of the 
        % enclosure! But, the main reason I am obssessed with the order 
        % here is that I want the scalar code to be consistent with the
        % vectorized version of this code as much as possible.
    end
    r(n-1,:) = (coeffs(n-1) - b(n-1,:)) + 2*(intval(t).'.*b(n,:));
    r(n,:) = coeffs(n) - b(n,:);
    
elseif strcmp(prec, 'quad')
    %% Use quadruple prec in residual calculation:
    % It does NOT matter a lot to make this part of the code vectorized,
    % mainly because quadruple precision is slow for long input vectors 
    % anyways and I do not intend to adverise it for long vectors...
    
    vec = b;
    r = intval(zeros(n,m));
    for i=n-2:-1:1
        a4 = vec(i:i+2,:);
        b4 = [-ones(size(t)), 2*t, -ones(size(t))]'; % multiplication by 2 is error-free
        if ( isintval(a4) | isintval(b4) )
            % quadruple precision ENCLOSURE for interval data
            %r(i) = Dot_int(a4', b4, 1, coeffs(i)); 
            
            r(i,:) = Dot_int(a4', b4, 1, coeffs(i));    
            
        else
            % quadruple precision ENCLOSURE for point data
                for j=1:m
                    r(i,j) = Dot_(1, coeffs(i), 2*t(j),b(i+1,j), ...
                    -1, b(i+2,j), -1, b(i,j), -2);
                end            
        end
    end
    if ( isintval(a4) | isintval(b4) )
        r(n-1) = Dot_int(vec(n-1:n)',[-1; 2*t], 1, coeffs(n-1));
        r(n) = Dot_int(vec(n)',-1, 1, coeffs(n));
    else
                r(n-1,j) = Dot_(1, coeffs(n-1), 2*t(j),b(n,j), ...
                    -1, b(n-1,j), -2);
                r(n,j) = Dot_(1, coeffs(n), -1, b(n,j), -2);
    end
end
end