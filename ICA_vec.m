function p= ICA_vec(t, intCoeffs)
% Naiive extension of Clenshaw to interval arithmetic. This is done e.g., 
% in Joldes' PhD thesis.
p = clenshaw_vec(intval(t), intCoeffs); 
end

function p = clenshaw_vec(x, coeffs)
% Clenshaw's algorithm to evaluate a polynomial p(x) at an INTERVAL x in 
% [-1, 1]. 
% x could be a vector.

bk1 = 0*x; % intervals
bk2 = bk1; % intervals
x2 = 2*x;  
len = size(coeffs,1); % degree of poly

for k = len:-1:1
    bk = coeffs(k) + x2.*bk1 - bk2;    
    bk2 = bk1; 
    bk1 = bk;
end
p = bk1 - x.*bk2;
end