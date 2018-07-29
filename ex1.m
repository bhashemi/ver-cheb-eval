%% Example 1, first part
%clc, 
clear all, close all, 
rng(1)
FS = 'fontsize'; fs = 15;
LW = 'LineWidth'; lw = 4;
loyolagray = 1/255*[200,200,200];
format short e

%% Generate some Chebyshev coefficients
%% Non-random case: the exponential function
n = 15; % chebLength = Number of Chebyshev coefficients to be used in Clenshaw
coeffs = 2* besseli((0:n-1)',1);
coeffs(1) = coeffs(1)/2;

% Another example with similar outputs
% f = chebfun(@(x) cos(exp(x)));
% coeffs = f.coeffs

%% Prepare coefficients of different radii
intCoeffs = coeffs;
maxrad_coeffs = max(rad(intCoeffs))

%% Adjust size of coefficients for the VERIFYFFT in barycentric method
%f = chebfun(ff);
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
%lenBary = size(baryCoeffs,1);
% max(abs(baryCoeffs(1:size(intCoeffs,1)) - intCoeffs)) = 0

%% Generate evaluation points
numpts = 10;
rng(1), t = 2*rand(numpts,1)-1;
[~,ind] = sort(t);
t = t(ind);

l = size(t,1);
maxrad_pts = max(rad(t))
check = all(in(t,infsup(-1,1)))

%% Call different methods and create plots for the paper
%% Direct cos-acos enclosures
fprintf('---- direct cos-acos enclosures (3.1) --\n')
tic
px_cos_acos = d_cos_acos(intCoeffs,t);
t_cos_acos = toc;
ulp_cos_acos = ulpcount(px_cos_acos);

%% divide & conquer
fprintf('---- divide & conquer (3.4): multiplication formulas --\n')
tic
px_divCon = d_div_con(intCoeffs,t);
t_divCon = toc;
ulp_divCon = ulpcount(px_divCon);

%% verifylss
fprintf('---- linear system solver (4.2) ----\n')
tic
for k = 1:size(t,1)
    pV(k) = LssErrBnd_Clenshaw_scl(t(k),intCoeffs);
end
t_LssErrBnd = toc;
ulp_LssErrBnd = ulpcount(pV);

%% ICA
fprintf('---- ICA  (4.3)--\n')
tic,
pICA = ICA_vec(t, intCoeffs);
t_ICA = toc;
ulp_ICA = ulpcount(pICA);

%%
fprintf('---- ICDC double prec (4.4) --\n')
prec = 'double';
tic,
pICDC = ICDC_vec(t, intCoeffs, prec);
t_ICDC = toc;
ulp_ICDC = ulpcount(pICDC);

%%
fprintf('---- ICDC quad prec (4.4) --\n')
prec = 'quad';
tic,
pICDC_Quad = ICDC_vec(t, intCoeffs, prec);
t_ICDC_Quad = toc;
ulp_ICDC_Quad = ulpcount(pICDC_Quad);

%% parallelepiped method 
fprintf('---- vectorized parallelepiped method (4.10) :error form already --\n')
[p_para_vec, tPara_vec] = ICA_para_err_vec(t, intCoeffs);
ulp_para_vec = ulpcount(p_para_vec);

%%
fprintf('---- vectorized Lohner QR method (4.10) --\n')
[p_LohnerQR_vec, tLohnerQR_vec] = ICA_QR_err_vec(t, intCoeffs);
ulp_LohnerQR_vec = ulpcount(p_LohnerQR_vec);

%% 
fprintf('---- vectorized ICA-eig-err (4.12) --\n')
[p_ICA_eig_err_vec, t_ICA_eig_err_vec] = ICA_eig_err_vec(t, intCoeffs);
ulp_ICA_eig_err  = ulpcount(p_ICA_eig_err_vec); % Not very narrow especially

%% verified 2nd barycentric evaluation
fprintf('---- 2nd barycentric representation (5.2) --\n')
lenBary = size(baryCoeffs,1);
tic
[x, w] = verchebpts(lenBary);
fvals = vercoeffs2vals(baryCoeffs);
px_bary = ver_bary(t, fvals, x, w);
t_bary = toc;
ulp_bary = ulpcount(px_bary);

%% Assemble all the results and create a few plots
clc
str = {'d-cos-acos', 'd-div-con', 'LssErrBnd', 'ICA', 'ICDC 2', 'ICDC 4   ',...
    'ICA-para-err', 'ICA-QR-err', 'ICA-eig-err',  'bary'};

radsNew = [rad(px_cos_acos)  rad(px_divCon')  rad(pV')  rad(pICA)  ...
    rad(pICDC')  rad(pICDC_Quad')...
    rad(p_para_vec')    rad(p_LohnerQR_vec') rad(p_ICA_eig_err_vec')   ...
    rad(px_bary)];

avgDigits = [mean(-log10(rad(px_cos_acos)))  mean(-log10(rad(px_divCon))) ...
    mean(-log10(rad(pV)))  mean(-log10(rad(pICA)))  ...
    mean(-log10(rad(pICDC)))   mean(-log10(rad(pICDC_Quad)))   ...
    mean(-log10(rad(p_para_vec)))  mean(-log10(rad(p_LohnerQR_vec)))  ...
    mean(-log10(rad(p_ICA_eig_err_vec)))   mean(-log10(rad(px_bary)))];

times = [t_cos_acos  t_divCon  t_LssErrBnd  t_ICA  t_ICDC t_ICDC_Quad ...
    tPara_vec  tLohnerQR_vec   t_ICA_eig_err_vec  t_bary];

format bank
accurate = avgDigits >= 15;

ulpsNew = [ulp_cos_acos  ulp_divCon  ulp_LssErrBnd  ulp_ICA,  ulp_ICDC ...
    ulp_ICDC_Quad  ulp_para_vec   ulp_LohnerQR_vec  ulp_ICA_eig_err   ulp_bary]'


subplot(1,3,1)
semilogy1 = semilogy(mid(t),radsNew,LW,lw);
set(semilogy1(1),'DisplayName','d-cos-acos','Marker','x','Color',[1 0 0],...%red
    'LineWidth',lw);

set(semilogy1(2),'DisplayName','d-div-con','Marker','<','LineStyle','-.',...
    'Color','m');
set(semilogy1(3),'DisplayName','LssErrBnd','color', 'blue');
set(semilogy1(4),'DisplayName','ICA','color', [0.3 0.7 0.9]);
set(semilogy1(5),'DisplayName','ICDC 2','Marker','o', 'color', 'c');
set(semilogy1(6),'DisplayName','ICDC 4','Marker','o', 'color', loyolagray)
set(semilogy1(7),'DisplayName','ICA-para-err','Marker','o','Color',[1 1 0],... % yellow
    'LineWidth',lw);
set(semilogy1(8),'DisplayName','ICA-QR-err','Marker','p','Color',[0.8 0.8 0.2],...
    'LineWidth',lw);
set(semilogy1(9),'DisplayName','ICA-eig-err', 'Color', [0 1 0]); % green
set(semilogy1(10),'DisplayName','bary','Marker','x','Color','k',...
    'LineWidth',lw);
xlabel('x', FS, fs)
ylabel('radius', FS, fs)

ax = subplot(1,3,2);
i = 1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'r');
%d-cos-acos

hold on
i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'm');
%d-div-con'

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor','b');
%LssErrBnd

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', [0.3 0.7 0.9]);
%ICA

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'c');
%ICDC 2

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor',loyolagray);
%ICDC 4

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'y');
%'para',

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', [0.8 0.8 0.2]);
%LohnerQR-vec

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor',[0 1 0]);
%ICA_eig_err

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'k');
%bary

ylim([0 16.1])
ylabel('average correct digits', FS, fs)
ax = gca;
set(ax, 'XTickLabel', {'d-cos-acos', 'd-div-con', 'LssErrBnd', 'ICA', ...
    'ICDC 2', 'ICDC 4', 'ICA-para-err', 'ICA-QR-err', ...
    'ICA-eig-err',  'bary'}, FS, fs-3)
set(gca,'xtick',1:i);
ax.XTickLabelRotation = 90;
xlim([0,i+1])
% Remove unnecessary parts of the second subplot
posDigInd = find( avgDigits .* (avgDigits>0) ); %index of avgDigits which are positive
m = min(avgDigits(posDigInd));
M = max(avgDigits(posDigInd));
ylim([m-0.1, M+0.1])

subplot(1,3,3)
i = 1;
h=bar(i,times(i));
set(h,'FaceColor', 'r');
%d-cos-acos

hold on
i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', 'm');
%d-div-con'

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor','b');
%LssErrBnd

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', [0.3 0.7 0.9]);
%ICA

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', 'c');
%ICDC 2

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor',loyolagray);
%ICDC 4

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', 'y');
%ICA-para-err

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', [0.8 0.8 0.2]);
%ICA-QR-err

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor',[0 1 0]);
%ICA_eig_err

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', 'k');
%bary

xlim([0,i+1])
ylabel('time (sec)', FS, fs)
ax = gca;
set(ax, 'XTickLabel', {'d-cos-acos', 'd-div-con', 'LssErrBnd', 'ICA', ...
    'ICDC 2', 'ICDC 4', 'ICA-para-err', 'ICA-QR-err', ...
    'ICA-eig-err', 'bary'}, FS, fs-3)

xTick=get(gca,'xtick');
xMax=max(xTick);
xMin=min(xTick);
newXTick=linspace(xMin,xMax,i);
set(gca,'xtick',1:i);
ax.XTickLabelRotation = 90;
text(-27,0.16,['degree = ' num2str(size(intCoeffs,1)-1) ',  rad(c) = ' num2str(maxrad_coeffs,2) ',    l = ' num2str(l) ',   rad(x) = ' num2str(maxrad_pts,2)], FS, fs)
shg
print(gcf,'-depsc','/Users/user/Desktop/My work/git/ver-cheb-eval/draft/figures/ex1');    