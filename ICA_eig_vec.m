function [pTrans, tTrans] = ICA_eig_vec(t, intCoeffs)
% Spectral transformation of interval Clenshaw algorithm.
% This is a vectorized code.

tic
j = size(t,1);
len = size(intCoeffs,1);

VI1 = [ 1; -1] * (-1i./(2*sqrt(1-intval(t).^2))).';
cTrans = kron(intCoeffs,VI1); % no for-loops here!

Dnew = [(t + 1i*sqrt(1-intval(t).^2)).'; (t - 1i*sqrt(1-intval(t).^2)).'];

btilde = zeros(2,j);
for k = len:-1:1
    btilde = Dnew .* btilde + cTrans(2*k-1:2*k,:);
end
bTrans = [(t + 1i*sqrt(1-intval(t).^2)).' .* btilde(1,:) + ...
          (t - 1i*sqrt(1-intval(t).^2)).' .* btilde(2,:);
          btilde(1,:) + btilde(2,:)];

pTrans = bTrans(1,:) - bTrans(2,:).*(t.');
if isreal(intCoeffs)
    pTrans = real(pTrans);
end
pTrans = pTrans.';
tTrans = toc;
end