function [p, tPara] = ICA_para_err_vec(t, intCoeffs)
% Vectorized version of the parallelepiped method applied to Clenshaw 
% algorithm as a discrete dynamical system.

tic
j = size(t,1);
M1 = 2*intval(t.');
Mid1 = mid(M1); 
midT = mid(t.');

bTilde = zeros(2,j); 
z = intval(zeros(2,j));

% Initial basis matrix
bOld1 = ones(1,j); bOld2 = zeros(1,j);
bOld3 = zeros(1,j); bOld4 = ones(1,j);
%Bold = [bOld1   bOld3;  bOld2   bOld4]

len = size(intCoeffs,1);
for k = len:-1:1
    cHat = [intCoeffs(k) .* ones(1,j); zeros(1,j)];
    
    %% Compute an approximation bTildeNew whose error we will enclose.
    bTildeNew = [bTilde(1,:) .* (2*midT);
               bTilde(1,:)] + bTilde(2,:) .* [-1; 0] + mid(cHat); 
                % floating point iteration for the approximation bTildeNew.
           
    %% Enclose the residual
    % Crucial to use 2*intval(t)= M1 in rHat, and NOT just 2*t:
    rHat = cHat + [bTilde(1,:) .* M1; bTilde(1,:)] + ...        
        intval(bTilde(2,:)) .* [-1; 0] - bTildeNew;  % Enclosure of 
    % residual in the approximation bTildeNew of (4.5). This should be a 
    % small quantity.
                                            
    %% Enclose the error
    % Choice of the basis matrix
    B1 = Mid1 .* bOld1 - bOld2;
    B2 = bOld1;
    B3 = Mid1 .* bOld3 - bOld4;
    B4 = bOld3;    
    %B = [B1   B3;  B2   B4]
    
    % Enclose inverse of B
    %B = [B1  B3; B2  B4];
    [BInv1, BInv2, BInv3, BInv4] = exactInv(B1, B2, B3, B4);
    
    % Form the parallelepiped
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
eHat(1,:) = B1 .* z(1,:) + B3 .* z(2,:);
eHat(2,:) = B2 .* z(1,:) + B4 .* z(2,:);

pApprox = bTildeNew(1,:) - bTildeNew(2,:) .* intval(t.'); % should consider
% errors in this additional computation as bTildeNew is not an interval.

p = pApprox + (eHat(1,:) - eHat(2,:) .* t.');
tPara = toc;
end

function [AI1, AI2, AI3, AI4] = exactInv(A1, A2, A3, A4)
% If A = [ A1, A3; A2, A4], then
% AInv = inv(A) = 1/(A1*A4 - A3*A2) * [A4, -A3; -A2,  A1].

%oneOverdetInv = 1./(A1 .* intval(A4) - A3 .* intval(A2));
oneOverdetInv = 1./(intval(A1) .* A4 - A3 .* intval(A2));
AI1 = oneOverdetInv .* A4;
AI2 = -oneOverdetInv .* A2;
AI3 = -oneOverdetInv .* A3;
AI4 = oneOverdetInv .* A1;
% AInv = [AI1   AI3; AI2   AI4]
end