%% Second part of Example 1 and Example 2
clc, 
clear all, close all, 
rng(1)
FS = 'fontsize'; fs = 15;
LW = 'LineWidth'; lw = 4;
loyolagray = 1/255*[200,200,200];
format short e

%% the exponential function
n = 15; % chebLength = Number of Chebyshev coefficients to be used in Clenshaw
coeffs = 2* besseli((0:n-1)',1);
coeffs(1) = coeffs(1)/2;

%% Generate coefficients of different radii
%intCoeffs = coeffs;
intCoeffs = midrad(coeffs, 100*eps);

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
numpts = 1000;
rng(1), t = 2*rand(numpts,1)-1; % bary is the best with this choice of points
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

%% divide & conquer
fprintf('---- divide & conquer (3.4): multiplication formulas --\n')
tic
px_divCon = d_div_con(intCoeffs,t);
t_divCon = toc;

%% ICA
fprintf('---- ICA  (4.3)--\n')
tic,
pICA = ICA_vec(t, intCoeffs);
t_ICA = toc;

%%
fprintf('---- ICDC double prec (4.4) --\n')
prec = 'double';
tic,
pICDC = ICDC_vec(t, intCoeffs, prec);
t_ICDC = toc;

%% parallelepiped method 
fprintf('---- vectorized parallelepiped method (4.10) :error form already --\n')
[p_para_vec, tPara_vec] = ICA_para_err_vec(t, intCoeffs);

%%
fprintf('---- vectorized Lohner QR method (4.10) --\n')
[p_LohnerQR_vec, tLohnerQR_vec] = ICA_QR_err_vec(t, intCoeffs);

%% 
fprintf('---- vectorized ICA-eig-err (4.12) --\n')
[p_ICA_eig_err_vec, t_ICA_eig_err_vec] = ICA_eig_err_vec(t, intCoeffs);

%% verified 2nd barycentric evaluation
fprintf('---- 2nd barycentric representation (5.2) --\n')
lenBary = size(baryCoeffs,1);
tic
[x, w] = verchebpts(lenBary);
fvals = vercoeffs2vals(baryCoeffs);
px_bary = ver_bary(t, fvals, x, w);
t_bary = toc;

%% Assemble all the results and create a few plots
str = {'d-cos-acos', 'd-div-con', 'ICA', 'ICDC 2',...
    'ICA-para-err', 'ICA-QR-err', 'ICA-eig-err',  'bary'};

radsNew = [rad(px_cos_acos)  rad(px_divCon')  rad(pICA)  ...
    rad(pICDC') rad(p_para_vec')    rad(p_LohnerQR_vec') ...
    rad(p_ICA_eig_err_vec')  rad(px_bary)];

avgDigits = [mean(-log10(rad(px_cos_acos)))  mean(-log10(rad(px_divCon))) ...
    mean(-log10(rad(pICA)))  mean(-log10(rad(pICDC)))  ...
    mean(-log10(rad(p_para_vec)))  mean(-log10(rad(p_LohnerQR_vec))) ...
    mean(-log10(rad(p_ICA_eig_err_vec)))  mean(-log10(rad(px_bary)))];

times = [t_cos_acos  t_divCon  t_ICA  t_ICDC  tPara_vec  tLohnerQR_vec ...
    t_ICA_eig_err_vec  t_bary];

close all
FigH = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
subplot(1,3,1)
semilogy1 = semilogy(mid(t),radsNew,LW,lw);
set(semilogy1(1),'DisplayName','cos-acos','Marker','x','Color',[1 0 0],...%red
    'LineWidth',lw);

set(semilogy1(2),'DisplayName','div-con','Marker','<','LineStyle','-.',...
    'Color','m');
set(semilogy1(3),'DisplayName','ICA','color', [0.3 0.7 0.9]);
set(semilogy1(4),'DisplayName','ICDC 2','Marker','o', 'color', 'c');
set(semilogy1(5),'DisplayName','ICA-para-err','Marker','o','Color',[1 1 0],... % yellow
    'LineWidth',lw);
set(semilogy1(6),'DisplayName','ICA-QR-err','Color',[0.8 0.8 0.2],...
    'LineWidth',lw);
set(semilogy1(7),'DisplayName','ICA-eig-err', 'Color', [0 1 0]); % green
set(semilogy1(8),'DisplayName','bary','Marker','x','Color','k',...
    'LineWidth',lw);
%set(semilogy1(7),'DisplayName','MVF', 'color', 'y');
xlabel('x', FS, fs)
ylabel('radius', FS, fs)
leg1 = legend(str, 'location', 'best');
set(leg1,...
    'Position',[0.1482143111527 0.692846347607066 0.149107142857143 0.182619647355164]);

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
set(h,'FaceColor', [0.3 0.7 0.9]);
%ICA

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'c');
%ICDC 2

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'y');
%ICA-para-err

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', [0.8 0.8 0.2]);
%ICA-QR-err-vec

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor',[0 1 0]);
%ICA-eig-err

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'k');
%bary

ylim([0 16.1])
ylabel('average correct digits', FS, fs)
ax = gca;
set(ax, 'XTickLabel', {'d-cos-acos', 'd-div-con', 'ICA', ...
    'ICDC 2', 'ICA-para-err', 'ICA-QR-err', ...
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
%cos-acos

hold on
i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', 'm');
%'div-con'

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
set(h,'FaceColor', 'y');
%'para',

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', [0.8 0.8 0.2]);
%LohnerQR-vec

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor',[0 1 0]);
%ICA-eig-err

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', 'k');
%bary

xlim([0,i+1])
ylabel('time (sec)', FS, fs)
ax = gca;
set(ax, 'XTickLabel', {'d-cos-acos', 'd-div-con', 'ICA', ...
    'ICDC 2', 'ICA-para-err', 'ICA-QR-err', ...
    'ICA-eig-err', 'bary'}, FS, fs-3)

xTick=get(gca,'xtick');
xMax=max(xTick);
xMin=min(xTick);
newXTick=linspace(xMin,xMax,i);
set(gca,'xtick',1:i);
ax.XTickLabelRotation = 90;
text(-22.5,0.19,['degree = ' num2str(size(intCoeffs,1)-1) ',  rad(c) = ' num2str(maxrad_coeffs,2) ',    l = ' num2str(l) ',   rad(x) = ' num2str(maxrad_pts,2)], FS, fs)
shg
print(gcf,'-depsc','/Users/user/Desktop/My work/git/ver-cheb-eval/draft/figures/ex346_new');