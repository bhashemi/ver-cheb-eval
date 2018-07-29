function values = vercoeffs2vals(coeffs)
%VERCOEFFS2VALS   Converts Chebyshev coefficients to enclosures for the 
% values at Chebyshev points of the 2nd kind.
% Equivalent to a verified DCT of type I.
% Adapted from Chebfun/chebtech2/vals2coeffs.

% Get the length of the input:
n = size(coeffs, 1);

% Trivial case (constant or empty):
if ( n <= 1 )
    values = coeffs; 
    return
end

% check for symmetry
isEven = max(abs(mid(coeffs(2:2:end,:))),[],1) == 0;
isOdd = max(abs(mid(coeffs(1:2:end,:))),[],1) == 0;

% Scale them by 1/2:
coeffs(2:n-1,:) = coeffs(2:n-1,:)/2; % error-free

% Mirror the coefficients (to fake a DCT using an FFT):
tmp = [coeffs; 
       coeffs(n-1:-1:2,:)];   
   
if ( isreal(coeffs) )
    % Real-valued case:
    values = real(verifyfft(tmp));
elseif ( isreal(1i*coeffs) )
    % Imaginary-valued case:
    values = 1i*real(verifyfft(imag(tmp)));
else
    % General case:
    values = verifyfft(tmp);
end

% Flip and truncate:
values = values(n:-1:1,:);

% enforce symmetry
values(:,isEven) = (values(:,isEven) + flipud(values(:,isEven)))/2;
values(:,isOdd) = (values(:,isOdd) - flipud(values(:,isOdd)))/2;

end