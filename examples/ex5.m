%% ex7.m
%clc, 
clear all, close all, 
rng(1)
FS = 'fontsize'; fs = 15;
LW = 'LineWidth'; lw = 4;
loyolagray = 1/255*[200,200,200];
format long e

%% Produce some Chebyshev coefficients
f = chebfun(@(x) sin(1./(x+1.003)));

%% Prepare coefficients of different radii
coeffs = f.coeffs;

intCoeffs = coeffs;
%intCoeffs = midrad(coeffs, 100*eps);
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

%% Prepare evaluation points
numpts = 10;
rng(1), t = 2*rand(numpts,1)-1;
[~,ind] = sort(t);
t = t(ind);

l = size(t,1); % length of the vector of evaluation points
maxrad_pts = max(rad(t))
check = all(in(t,infsup(-1,1)))
output_rad_expect = len * max(rad(t)) % ? max(rad(intCoeffs))?

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

%% parallelepiped method 
fprintf('---- vectorized parallelepiped method (4.10) :error form already --\n')
[p_para_vec, tPara_vec] = ICA_para_err_vec(t, intCoeffs);

%%
fprintf('---- vectorized Lohner QR method (4.10) --\n')
[p_LohnerQR_vec, tLohnerQR_vec] = ICA_QR_err_vec(t, intCoeffs);

%% Vectorized error form of TICA
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

%% verifylss
fprintf('---- linear system solver ----\n')
tic
for k = 1:size(t,1)
    pV(k) = LssErrBnd_Clenshaw_scl(t(k),intCoeffs);
end
t_LssErrBnd = toc;

%% Assemble all the results and create a few plots
str = {'d-cos-acos', 'd-div-con', 'LssErrBnd', 'ICA-para-err', 'ICA-QR-err', ...
    'ICA-eig-err', 'bary'};

radsNew = [rad(px_cos_acos)  rad(px_divCon')  rad(pV)' rad(p_para_vec')...
    rad(p_LohnerQR_vec') rad(p_ICA_eig_err_vec')  rad(px_bary)];

avgDigits = [mean(-log10(rad(px_cos_acos))) mean(-log10(rad(px_divCon)))...
    mean(-log10(rad(pV)))  mean(-log10(rad(p_para_vec)))  ...
    mean(-log10(rad(p_LohnerQR_vec))) mean(-log10(rad(p_ICA_eig_err_vec))) ...
     mean(-log10(rad(px_bary)))];

format long e
times = [t_cos_acos   t_divCon   t_LssErrBnd   tPara_vec  ...
    tLohnerQR_vec   t_ICA_eig_err_vec   t_bary];

close all
FigH = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
subplot(1,3,1)
semilogy1 = semilogy(mid(t),radsNew,LW,lw);
set(semilogy1(1),'DisplayName','d-cos-acos','Marker','x','Color',[1 0 0],...%red
    'LineWidth',lw);
set(semilogy1(2),'DisplayName','d-div-con','Marker','<','LineStyle','-.',...
    'Color','m');
set(semilogy1(3),'DisplayName','LssErrBnd','color', 'blue');

set(semilogy1(4),'DisplayName','ICA-para-err','Marker','o','Color',[1 1 0],... % yellow
    'LineWidth',lw);

set(semilogy1(5),'DisplayName','ICA-QR-err','Marker','p','Color',[0.8 0.8 0.2],...
    'LineWidth',lw);
set(semilogy1(6),'DisplayName','ICA-eig-err', 'Color', [0 1 0]); % green
set(semilogy1(7),'DisplayName','bary','Marker','x','Color','k',...
    'LineWidth',lw);

xlabel('x', FS, fs)
ylabel('radius', FS, fs)
leg1 = legend(str, 'location', 'ne');
set(leg1,...
    'Position',[0.1482143111527 0.712846347607066 0.149107142857143 0.182619647355164]);


ax = subplot(1,3,2);
i = 1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'r');
%cos-acos

hold on
i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'm');
%'div-con'

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor','b');

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
%ICA-eig-err


i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', 'k');
%bary

ylim([0 16.1])
ylabel('average correct digits', FS, fs)
ax = gca;
set(ax, 'XTickLabel', {'d-cos-acos', 'd-div-con', 'LssErrBnd',   ...
    'ICA-para-err', 'ICA-QR-err', 'ICA-eig-err',  'bary'}, FS, fs-3)

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
set(h,'FaceColor','b');

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
ylim([0, max(times)+0.3])
ylabel('time (sec)', FS, fs)
ax = gca;
xticklabels(ax, {'d-cos-acos', 'd-div-con', 'LssErrBnd', 'ICA-para-err', ...
    'ICA-QR-err', 'ICA-eig-err',  'bary'});

xTick=get(gca,'xtick'); 
xMax=max(xTick);
xMin=min(xTick);
newXTick=linspace(xMin,xMax,i);
set(gca,'xtick',1:i);
ax.XTickLabelRotation = 90;
text(-19,49.5,['degree = ' num2str(size(intCoeffs,1)-1) ',  rad(c) = ' num2str(maxrad_coeffs,2) ',    l = ' num2str(l) ',   rad(x) = ' num2str(maxrad_pts,2)], FS, fs)
