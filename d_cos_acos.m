function px = d_cos_acos(coeffs, x)
% Compute enclosures for values of the Chebyshev expansion 
% p(x) = coeffs(1) T_0(x) + coeffs(1) T_1(x) + ... + coeffs(n) T_{n-1}(x),
% where coeffs is an n x 1 column vector of Chebyshev coefficients and 
% x is an m x 1 column vector of evaluation points.

n = size(coeffs,1);
k = 0:n-1;  % row vector
Tx = cos(k.*acos(intval(x))); % m x n matrix
px = Tx*coeffs;               % m x 1 column vector of output enclosures
end