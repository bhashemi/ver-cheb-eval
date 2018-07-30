%% Forward error analysis of Chebfun/feval
clc, clear all, close all, 
rng(1)
FS = 'fontsize'; fs = 15;
LW = 'LineWidth'; lw = 4;
MS = 'markersize'; ms = 10;
loyolagray = 1/255*[200,200,200];
format short e
dom = [-1, 1];

%% Generate some Chebyshev coefficients
f = cheb.gallery('sinefun1')
%dom = f.domain
coeffs = f.coeffs;
intCoeffs = coeffs;
maxrad_coeffs = max(rad(intCoeffs))

%% Adjust size of coefficients for the VERIFYFFT in barycentric method
len = size(intCoeffs,1)
%coeffs = [f.coeffs; zeros(2^nextpow2(len)-len,1)];
% VERCOEFFS2VALS (and COEFFS2VALS) first converts an input vector of size m
% to an vector of size 2m-2. 
% COEFFS2VALS then calls FFT. Analogously, VERCOEFFS2VALS calls INTLAB's
% VERIFYFFT on the vector of size 2m-2. But, VERIFYFFT (in INTLAB Version 
% 10) works only if 2m-2 is a power of two. To keep things simple, we first
% prolong the coefficients of the Chebfun object so that 2*(its length) - 2
% is a power of two!
newSize = (2^nextpow2(2*len-2) + 2)/2;
dif = newSize - len;
baryCoeffs = [intCoeffs; zeros(dif,1)];

%% Choose evaluation points
numpts = 100;
t = 2*rand(numpts,1)-1;
[~,ind] = sort(t);
t = t(ind);

%% map the points to [-1,1]
if dom(1) ~= -1 && dom(2) ~= 1
    dom = intval(dom);
   t = (2*t - (dom(2) + dom(1)))./(dom(2) - dom(1));
end

%%
l = size(t,1); % length of the vector of evaluation points
maxrad_pts = max(rad(t))
check = all(in(t,infsup(-1,1)))

chebvals = f(t);

%% Call different methods and create plots for the paper
fprintf('---- vectorized affine transformation --\n')
[p_ICA_eig, t_ICA_eig] = ICA_eig_vec(t, intCoeffs);

ins = in(chebvals,p_ICA_eig);
sum_ins_cheb_trans = sum(ins)

%% vectorized ICA
fprintf('---- vectorized ICA --\n')
tic,
pICA = ICA_vec(t, intCoeffs);
t_ICA = toc;

%% Direct cos-acos enclosures
fprintf('---- direct cos-acos enclosures --\n')
tic
px_cos_acos = d_cos_acos(intCoeffs,t);
t_cos_acos = toc;

%% verified barycentric evaluation
fprintf('---- barycentric representation --\n')
lenBary = size(baryCoeffs,1);
tic
[x, w] = verchebpts(lenBary);
fvals = vercoeffs2vals(baryCoeffs);
px_bary = ver_bary(t, fvals, x, w);
t_bary = toc;

%% divide & conquer
fprintf('---- divide & conquer: multiplication formulas --\n')
tic
px_divCon = d_div_con(intCoeffs,t);
t_divCon = toc;

in_bary_cosAcos = sum(in(px_bary,px_cos_acos))
in_bary_divCon = sum(in(px_bary,px_divCon'))
format long

%% Intersection of all could be adviceable? a heurisitic... no method is 
% always great!
tic
px_intersect = intersect(px_bary, px_cos_acos);
t_px_intersect = toc + t_bary + t_cos_acos;

%% vectorized ICDC double
fprintf('---- ICDC double prec --\n')
prec = 'double';
tic,
pICDC = ICDC_vec(t, intCoeffs, prec);
t_ICDC = toc;
ins_double = in(chebvals,pICDC');
sum_ins_cheb_double = sum(ins_double)

%%
fprintf('---- ICDC quadruple prec --\n')
prec = 'quad';
tic,
pICDCQuad = ICDC_vec(t, intCoeffs, prec);
t_ICDCQuad = toc;
ins_Quad = in(chebvals,pICDCQuad');
sum_ins_cheb_Quad = sum(ins_Quad); 

t(60)
chebvals(60)
pICDCQuad(60)

% Display more digitis of the results
c = vpa(chebvals(60),20)
% c =
% 2.4456137878210535419
x = vpa(inf(pICDCQuad(60)),20)
% x =
% 2.4456137878210544301
y = vpa(sup(pICDCQuad(60)),20)
% y =
% 2.4456137878210548742

%%
plot(f)
hold on
plot(t,f(t), 'kx')
ind=find(~ins_Quad);
plot(t(ind), f(t(ind)), 'ro', MS, ms)
xlabel('x')
ylabel('y = f(x)')
set(gca, FS, 15)
leg = legend('sinefun1', 'floating point values', 'values outside enclosures', ...
    'orientation', 'horizontal', 'location', 'n');
set(leg, FS, 12.6)
set(leg,...
    'Position',[0.146428571428584 0.8738096015794 0.741071428571429 0.0464285714285715])
