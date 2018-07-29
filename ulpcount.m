function ulpB = ulpcount(B)
% Count an upper bound on the accuracy of enclosures, i.e., the number of
% representable machine numbers in the interior of the interval B.
% See Page 47 of Hammer, Hocks, Kulisch, Ratz, Numerical Toolbox for
% Verified Computing I, 1993.

[m, n] = size(B);
counts = zeros(m,n);
%dist = zeros(m,n);
for i=1:m
    for j=1:n
        % The matrix is symmetric, so only do the job for its upper
        % triangular part.
        infZ = inf(B(i,j));
        %getbits(infZ)
        ulp = 0;
        while infZ + eps(infZ) < sup(B(i,j))
            %binZ = binstr(infZ);  % binary representation of different infZ values 
            % could be displayed exactly the same; display their binary representations 
            % to see the difference.
            infZ = infZ + eps(infZ);
            %getbits(infZ)
            ulp = ulp+1;
        end
        counts(i,j) = ulp+1;
        %fprintf('\n The enclosure is accurate to %d ulps\n', counts(i,j));
    end    
end
ulpB = mean(counts(:));

end