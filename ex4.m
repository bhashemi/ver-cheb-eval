%% All three of the Examples 3, 4 and 5
clc, 
clear all, close all, 
rng(1)
FS = 'fontsize'; fs = 15;
LW = 'LineWidth'; lw = 4;
loyolagray = 1/255*[200,200,200];
format short e

%% Prepare coefficients of different radii
rng(1), f = randnfun(0.0069) % degree 999
coeffs = f.coeffs;
%intCoeffs = real(coeffs);
intCoeffs = midrad(real(coeffs), 100*eps);

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

%% Prepare evaluation points
numpts = 5e4;

rng(1), t = 2*rand(numpts,1)-1;
[~,ind] = sort(t);
t = t(ind);

t = midrad(mid(t), 1e2*eps);

l = size(t,1);
maxrad_pts = max(rad(t))
check = all(in(t,infsup(-1,1)))

%% Call different methods and create plots for the paper
%% verified barycentric evaluation
fprintf('---- barycentric representation --\n')
lenBary = size(baryCoeffs,1);
tic
[x, w] = verchebpts(lenBary);
fvals = vercoeffs2vals(baryCoeffs);
px_bary = ver_bary(t, fvals, x, w);
t_bary = toc;

%% Direct cos-acos enclosures
fprintf('---- direct cos-acos enclosures --\n')
tic
px_cos_acos = d_cos_acos(intCoeffs, t);
t_cos_acos = toc;

%% divide & conquer
fprintf('---- divide & conquer: multiplication formulas --\n')
tic
px_divCon = d_div_con(intCoeffs,t);
t_divCon = toc;

in_bary_cosAcos = sum(in(px_bary,px_cos_acos))
in_bary_divCon = sum(in(px_bary,px_divCon'))

%% 
fprintf('---- Vectorized error form of TICA 2--\n')
[p_ICA_eig_err_vec, t_ICA_eig_err_vec] = ICA_eig_err_vec(t, intCoeffs);

%% parallelepiped method 
fprintf('---- Vectorized parallelepiped method (error form already) --\n')
[p_para_vec, t_para_vec] = ICA_para_err_vec(t, intCoeffs);

%%
fprintf('---- Lohner''s  QR method, vectorized --\n')
tic
[p_LohnerQR_vec, tLohnerQR_vec] = ICA_QR_err_vec(t, intCoeffs);

%%
str = {'d-cos-acos', 'd-div-con', ...
    'ICA-para-err', 'ICA-QR-err', 'ICA-eig-err',...
    'bary'};
rads = [rad(px_cos_acos)  rad(px_divCon)'  rad(p_para_vec') ...
    rad(p_LohnerQR_vec')  rad(p_ICA_eig_err_vec')  rad(px_bary)];

avgDigits = [mean(-log10(rad(px_cos_acos))) mean(-log10(rad(px_divCon))) ...
    mean(-log10(rad(p_para_vec)))  mean(-log10(rad(p_LohnerQR_vec)))  ...
    mean(-log10(rad(p_ICA_eig_err_vec))) mean(-log10(rad(px_bary))) ];

times = [t_cos_acos  t_divCon  t_para_vec      tLohnerQR_vec  ...
    t_ICA_eig_err_vec  t_bary];

close all
FigH = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
subplot(1,3,1)

h1 = semilogy(mid(t),rads(:,1), 'Marker','x', 'color', 'r',LW,lw);
hold on
h2 = semilogy(mid(t),rads(:,2), 'm<-.',LW,lw);
h3 = semilogy(mid(t),rads(:,3), 'o-y',LW,lw);
h4 = semilogy(mid(t),rads(:,4), 'Marker','o', 'color', [0.8 0.8 0.2], LW,lw);
h5 = semilogy(mid(t),rads(:,5), 'Marker','p', 'color', 'g',LW,lw);
h6 = semilogy(mid(t),rads(:,6), 'Marker','x', 'color', 'k',LW,lw);
h1 = semilogy(mid(t),rads(:,1), 'Marker','x', 'color', 'r',LW,lw);
h2 = semilogy(mid(t),rads(:,2), 'm<-.',LW,lw);
h3 = semilogy(mid(t),rads(:,3), 'o-y',LW,lw);
h4 = semilogy(mid(t),rads(:,4), 'Marker','o', 'color', [0.8 0.8 0.2], LW,lw);
h5 = semilogy(mid(t),rads(:,5), 'Marker','p', 'color', 'g',LW,lw);

xlabel('x', FS, fs)
ylabel('radius', FS, fs)

%plot ICA_eig_err again as it is below bary:
legend([h1 h2 h3 h4 h5 h6], str)
leg1 = legend(str, 'location', 'best');
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
ax = gca;
set(ax, 'XTickLabel', {'d-cos-acos', 'd-div-con', 'ICA-para-err', 'ICA-QR-err', ...
    'ICA-eig-err',  'bary'}, FS, fs-3)
xlim([0,i+1])
ylabel('average correct digits', FS, fs)
% Remove unnecessary parts of the second subplot
posDigInd = find( avgDigits .* (avgDigits>0) ); %index of avgDigits which are positive
m = min(avgDigits(posDigInd));
M = max(avgDigits(posDigInd));
ylim([m-0.1, M+0.1])
set(ax,'xtick',1:i);
ax.XTickLabelRotation = 90;


ax = subplot(1,3,3);
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
set(h,'FaceColor', 'y');
%ICA-para-err

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', [0.8 0.8 0.2]);
%ICA-QR-err-vec

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor',[0 1 0]);
%ICA_eig_err

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', 'k');
%bary

ax = gca;
set(ax, 'XTickLabel', {'d-cos-acos', 'd-div-con', 'ICA-para-err', 'ICA-QR-err', ...
    'ICA-eig-err',  'bary'}, FS, fs-3)
ax.XTickLabelRotation = 90;
xlim([0,i+1])
ylabel('time (sec)', FS, fs)

set(ax,'xtick',1:i);
ax.XTickLabelRotation = 90;

text(-18,0.85,['degree = ' num2str(size(intCoeffs,1)-1) ',  rad(c) = ' num2str(maxrad_coeffs,2) ',    l = ' num2str(l) ',   rad(x) = ' num2str(maxrad_pts,2)], FS, fs)
shg
print(gcf,'-depsc','/Users/user/Desktop/My work/git/ver-cheb-eval/draft/figures/exdeg100_new');