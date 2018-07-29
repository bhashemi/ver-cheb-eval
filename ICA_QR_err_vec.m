function [p, tLohnerQR_vec] = ICA_QR_err_vec(x, intCoeffs)
% Vectorized version of the Lohner's QR method applied to Clenshaw 
% algorithm as a discrete dynamical system. 
% The input x is a column vector.

% We explicitly compute an interval matrix containing the exact QR 
% decomposition of 2 x 2 matrices (explicit in the sense that we don't call
% MATLAB's qr command as it does not make sense for vectors of evaluation 
% points). We then choose B (which is needed only approximately) to be the 
% midpoint of that interval matrix. We then need an enclosure for the exact
% inverse of mid(B), which simply belongs to B.'. So, there is no need to
% matrix inversion or calling verifylss.

tic
M1 = 2*intval(x.'); % Entry (1,1) of M.
j = size(x,1);

Mid1 = mid(M1);
%midA = [Mid1   Mid3;  
%        Mid2   Mid4]
midX = mid(x.');

bTilde = zeros(2,j); z = intval(zeros(2,j)); % Error in \tilde{y}_0 which is zero because \tilde{y}_0 is exactly representable.

% Initial basis matrix: identity 
bOld1 = ones(1,j); bOld2 = zeros(1,j);
bOld3 = zeros(1,j); bOld4 = ones(1,j);

len = size(intCoeffs,1);
for k = len:-1:1
    cHat = [intCoeffs(k) .* ones(1,j); zeros(1,j)];
    
    %% Compute an approximation bTildeNew whose error we will enclose.
    bTildeNew = [bTilde(1,:) .* (2*midX);
               bTilde(1,:)] + bTilde(2,:) .* [-1; 0] + mid(cHat);     
    
    % Crucial to use 2*intval(t)= M1 in rHat, and NOT just 2*t:
    rHat = cHat + [bTilde(1,:) .* M1; bTilde(1,:)] + ...        
        intval(bTilde(2,:)) .* [-1; 0] - bTildeNew;  % Enclosure of 
    % residual in the approximation bTildeNew of (4.5). This should be a 
    % small quantity.
                                            
    %% Choice of the basis matrix without calling MATLAB's qr command
    % We compute the exact B like a symbolic QR factorization of any 2x2
    % matrix so that BInv = B'
    %% Enclose the error
    % Choice of the basis matrix: BB is needed only    
    %BB = mid(A)*Bold;        % Store BB = [BB1   BB3;  BB2   BB4] to use 
    % Lohner's QR 
    BB1 = Mid1 .* bOld1 - bOld2;
    BB2 = bOld1;
    
    %% Compute an enclosure for the Q factor of the QR decomposition of BB    
    if in(zeros(2,j), intval([BB1; BB2]))
        % Q = I_{2 x 2}
        B_exact1 = intval(ones(1,j));
        B_exact2 = intval(zeros(1,j));
        B_exact3 = B_exact2;
        B_exact4 = B_exact1;
    else
        oneOverNrmCol1 = 1 ./ sqrt(intval(BB1).^2 + intval(BB2).^2);
        B_exact1 = oneOverNrmCol1 .* BB1;
        B_exact2 = oneOverNrmCol1 .* BB2;
        B_exact3 = -B_exact2; % B_exact3 = -oneOverNrmCol1 .* BB2;
        B_exact4 =  B_exact1; % B_exact4 =  oneOverNrmCol1 .* BB1;
    end

    %B = mid(B_exact);
    B1 = mid(B_exact1);
    B2 = mid(B_exact2); 
    B3 = mid(B_exact3); 
    B4 = mid(B_exact4);
    
    %% Enclose inverse of B. Don't need verifylss: B^{-1} \in [B_exact]^T.
    %BInv = B_exact.';
    BInv1 = B_exact1;
    BInv2 = B_exact3;
    BInv3 = B_exact2;
    BInv4 = B_exact4;
    
    %% Form the parallelepiped
    % Step 1: Compute T := M*Bold, where T = [T1  T3; T2  T4]
    T1 = M1 .* bOld1 - bOld2;
    T2 = bOld1;
    T3 = M1 .* bOld3 - bOld4;
    T4 = bOld3;
    
    % Step 2: Compute S := BInv * T, where S = [S1  S3; S2  S4]
    S1 = BInv1 .* T1 + BInv3 .* T2;
    S2 = BInv2 .* T1 + BInv4 .* T2;
    S3 = BInv1 .* T3 + BInv3 .* T4;
    S4 = BInv2 .* T3 + BInv4 .* T4;
    
    % Step 3: Compute s := S * z, where s = [s1; s2]
    s1 = S1 .* z(1,:) + S3 .* z(2,:);
    s2 = S2 .* z(1,:) + S4 .* z(2,:);
    
    % Step 4: Compute ss := BInv * rHat, where ss = [ss1; ss2]
    ss1 = BInv1 .* rHat(1,:) + BInv3 .* rHat(2,:);
    ss2 = BInv2 .* rHat(1,:) + BInv4 .* rHat(2,:);

    % Step 5: Update z := s + ss
    z(1,:) = s1 + ss1;
    z(2,:) = s2 + ss2;
    
    %% prepare for the next iterate
    bOld1 = B1; bOld2 = B2; bOld3 = B3; bOld4 = B4;
    bTilde = bTildeNew;
end
% Final enclosure for the error
eHat(1,:) = B1 .* z(1,:) + B3 .* z(2,:);
eHat(2,:) = B2 .* z(1,:) + B4 .* z(2,:);

pApprox = bTildeNew(1,:) - bTildeNew(2,:) .* intval(x.'); % should consider
% errors in this additional computation as bTildeNew is not an interval.

p = pApprox + (eHat(1,:) - eHat(2,:) .* x.');
tLohnerQR_vec = toc;
end