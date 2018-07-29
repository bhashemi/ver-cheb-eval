function px = d_div_con(coeffs, x)
% Divide and conquer idea to enclose Chebyshev polynomials using product
% formulas.
n = size(coeffs,1);

%% The following version is faster than the above.
k = floor(log2(n));
Tx(1,:) = intval(ones(size(x)));
Tx(2,:) = intval(x);
for i = 1:k-1 % We do the last time separately afterwards
    % N.B.: The for-loop is called only log_2(n) times!
    %Tx = div_con(Tx,x);
    Tx = div_con(Tx,x,2^i);
end

% ... afterwards, so that we do not compute unnecessary values from (2^k)+1
% all the way up to 2^(k+1). We just need to go for indices from (2^k)+1 to
% n.
num = numel(2^k+1:2:n);
num2 = numel(2^k+2:2:n);
Tx(2^k+1:2:n,:) = 2*Tx(2^(k-1)+1:1:2^(k-1)+ num,:).^2-1;
Tx(2^k+2:2:n,:) = 2*Tx(2^(k-1)+2:1:2^(k-1)+ 1 + num2,:) .* Tx(2^(k-1)+1:1:2^(k-1)+ num2,:) - x.';
px = coeffs'*Tx;
end

function Tx = div_con(Tx, x, k)
Tx(k+1:2:2*k,:) = 2*Tx(k/2+1:1:k,:).^2-1;
Tx(k+2:2:2*k,:) = 2*Tx(k/2+2:1:k+1,:) .* Tx(k/2+1:1:k,:) - x.';
end