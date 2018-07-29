function [p, b] = LssErrBnd_Clenshaw_scl(t,coeffs)
% %% Build the linear system and solve it in interval arithmetic
n = size(coeffs, 1);
col = intval(zeros(n,1)); col(1) = 1;
row = intval(zeros(1,n));
row(1) = 1; row(2) = -2*t; row(3) = 1;
A = toeplitz(col, row);
b = verifylss(A, coeffs);
p = b(1) - b(2)*t;
end