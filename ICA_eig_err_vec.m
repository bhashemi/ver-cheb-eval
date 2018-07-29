function [pTrans, tTrans] = ICA_eig_err_vec(t, intCoeffs)
% This is the error form of spectral transformation of interval 
% clenshaw algorithm".
%
% \hat{b_k} = M \hat{b_{k+1}} + \hat{c_k},    Formula (4.9) in the draft
%
% where a VECTOR of evaluation points t is given.
% The error form is Formula (4.10) in the draft. Specifically, here we
% compute an enclosure FOR THE ERROR in each term of the transformed 
% Clenshaw recurrence, rather than enclosing each term directly.

tic
t = intval(t);
midT = mid(t);
j = size(t,1);

bTilde = zeros(2,j);
eOld = zeros(2,j);
len = size(intCoeffs,1);

VI1 = [ 1; -1] * (-1i./(2*sqrt(1-intval(t).^2))).';
VI2 = [ ((1i*sqrt(1-intval(t).^2) - t) .* (-1i./(2*sqrt(1-intval(t).^2)))).'; 
    ((1i*sqrt(1-intval(t).^2) + t) .* (-1i./(2*sqrt(1-intval(t).^2)))).'];

Dnew = [(t + 1i*sqrt(1-intval(t).^2)).'; (t - 1i*sqrt(1-intval(t).^2)).'];

for k = len:-1:1
    bTildeNew = [mid(intCoeffs(k)) + (bTilde(1,:) .* (2*midT.') - bTilde(2,:));
               bTilde(1,:)];
    rHat = [intCoeffs(k) + (bTilde(1,:) .* (2*t.') - bTilde(2,:)) - bTildeNew(1,:); 
            bTilde(1,:) - intval(bTildeNew(2,:)) ];        
                                            
    rCheck = VI1 .* rHat(1,:) + VI2 .* rHat(2,:);

    eNew = Dnew .* eOld + rCheck;  % interval iteration for the error
                                                  
    %% Prepare for the next iterate:
    bTilde = bTildeNew;
    eOld = eNew;
end

eNew = [(t + 1i*sqrt(1-intval(t).^2)).' .* eNew(1,:) + ...
    (t - 1i*sqrt(1-intval(t).^2)).' .* eNew(2,:);
    eNew(1,:) + eNew(2,:)];

err = eNew(1,:) - eNew(2,:) .* (t.');

pTilde = bTildeNew(1,:) - bTildeNew(2,:) .* (t.');
pTrans = pTilde + err;
if isreal(intCoeffs)
    pTrans = real(pTrans);
end
tTrans = toc;

end