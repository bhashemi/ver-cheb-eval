function s = Dot_int(a,b, pm, c)
% Quadruple precision ENCLOSURE for a'*b pm c where a and b are interval 
% vectors and c is a scalar interval. This is done via case-distinction 
% S. M. Rump, Private communication, 2 Nov. 2017.

a = intval(a); 
b = intval(b); 
c = intval(c);

bb = b';
% case distinction as in Table 1 in p. 24 of Perspectives on Enclosure
% Methods.
con1467 = ( (inf(a) >= 0 & inf(bb) >= 0) |... % Conditions 1, 4, 6 and 7 in
    (sup(a) <= 0 & inf(bb) >= 0) | ...        % the table
    (sup(a) <= 0 & in0(0,bb) ) | ...
    (in0(0,a) & inf(bb) >= 0) );
aLbounds = con1467.*inf(a) + (~con1467).*sup(a);

con4567 = ( (sup(a) <= 0 & inf(bb) >= 0) | ... % Conditions 4, 5, 6 and 7 
    (sup(a) <= 0 & sup(bb) <=0 ) | ...         % in the table
    (sup(a) <= 0 & in0(0,bb) ) | ...
    (in0(0,a) & inf(bb) >= 0) );
bLbounds = con4567.*sup(bb) + (~con4567).*inf(bb);

con1347 = ( (inf(a) >= 0 & inf(bb) >= 0) |... % 1
    (inf(a) >= 0 & in0(0,bb) ) | ...          % 3
    (sup(a) <= 0 & sup(bb) <=0 ) | ...        % 4
    (in0(0,a) & inf(bb) >= 0) );              % 7
aUbounds = con1347.*sup(a) + (~con1347).*inf(a);

con1237 = ( (inf(a) >= 0 & inf(bb) >= 0) |... % 1
    (inf(a) >= 0 & sup(bb) <= 0 ) | ...       % 2
    (inf(a) >= 0 & in0(0,bb) ) | ...          % 3
    (in0(0,a) & inf(bb) >= 0) );              % 7
bUbounds = con1237.*sup(bb) + (~con1237).*inf(bb);
bLbounds = bLbounds';
bUbounds = bUbounds';

if ( any(in0(0,a)) & any(in0(0,b)) )
    error('Dot_int','should include the 9th case, too!')
end

if pm == 1
    sInf = Dot_(aLbounds,bLbounds, pm, inf(c), -2);
    sSup = Dot_(aUbounds,bUbounds, pm, sup(c), -2);
        
elseif pm == -1
    sInf = Dot_(aLbounds,bLbounds, pm, sup(c), -2);
    sSup = Dot_(aUbounds,bUbounds, pm, inf(c), -2);
end    
s = infsup(inf(sInf), sup(sSup));
end