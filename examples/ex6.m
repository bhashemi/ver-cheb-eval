%% ex9New.m
%clc, 
clear all, close all, 
rng(1)
FS = 'fontsize'; fs = 15;
LW = 'LineWidth'; lw = 4;
loyolagray = 1/255*[200,200,200];
format short e

%% Generate some Chebyshev coefficients
f = randnfun(0.0007)
coeffs = real(f.coeffs);
%% Prepare coefficients of different radii
%intCoeffs = coeffs;
intCoeffs = midrad(coeffs, 10*eps);
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

%% Prepare evaluation points
numpts = 1000;
rng(1), t = 2*rand(numpts,1)-1; 
[~,ind] = sort(t);
t = t(ind);
t = midrad(mid(t), 10*eps);
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

%% parallelepiped method 
% fprintf('---- vectorized parallelepiped method (4.10) :error form already --\n')
% [p_para_vec, tPara_vec] = ICA_para_err_vec(t, intCoeffs);

%%
% fprintf('---- vectorized Lohner QR method (4.10) --\n')
% [p_LohnerQR_vec, tLohnerQR_vec] = ICA_QR_err_vec(t, intCoeffs);

%%
fprintf('---- vectorized ICA-eig (4.11) --\n')
[p_ICA_eig_vec, t_ICA_eig_vec] = ICA_eig_vec(t, intCoeffs);
%Almost always worse than ICA_eig_err_vec!

%%
fprintf('---- vectorized ICA-eig-err (4.12) --\n')
%[p_ICA_eig_err_vec, t_ICA_eig_err_vec] = ICA_eig_err_col_vec(t, intCoeffs);
[p_ICA_eig_err_vec, t_ICA_eig_err_vec] = ICA_eig_err_vec(t, intCoeffs);

%% verified 2nd barycentric evaluation
fprintf('---- 2nd barycentric representation (5.2) --\n')
lenBary = size(baryCoeffs,1);
tic
[x, w] = verchebpts(lenBary);
fvals = vercoeffs2vals(baryCoeffs);
px_bary = ver_bary(t, fvals, x, w);
t_bary = toc;

%% Intersection of all could be adviceable? a heurisitic... no method is 
% always great!
tic
px_intersect = intersect(px_bary, px_cos_acos);
t_px_intersect = toc + t_bary + t_cos_acos;

%%
str1 = {'d-cos-acos', 'd-div-con', 'ICA-eig', 'ICA-eig-err', 'bary'};
str2 = {'d-cos-acos', 'd-div-con',  'ICA-eig', 'ICA-eig-err', 'bary', 'bary \cap d-cos-acos'};

rads = [rad(px_cos_acos)   rad(px_divCon)' rad(p_ICA_eig_vec) ...
    rad(p_ICA_eig_err_vec') rad(px_bary)   rad(px_intersect)];

avgDigits = [mean(-log10(rad(px_cos_acos))) ...
    mean(-log10(rad(px_divCon'))) mean(-log10(rad(p_ICA_eig_vec))) ...
    mean(-log10(rad(p_ICA_eig_err_vec))) mean(-log10(rad(px_bary)))...
    mean(-log10(rad(px_intersect)))];

times = [t_cos_acos  t_divCon    t_ICA_eig_vec  t_ICA_eig_err_vec ...
    t_bary  t_px_intersect];

%%
close all
%FigH = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
subplot(1,2,1)
h1 = semilogy(mid(t),rads(:,1), 'Marker','x', 'color', 'r',LW,lw);
hold on
h2 = semilogy(mid(t),rads(:,2), 'm<-.',LW,lw);
h3 = semilogy(mid(t),rads(:,3), 'color', [0.2 0.8 0] ,LW,lw);
h4 = semilogy(mid(t),rads(:,4), 'color', [0 1 0], LW,lw);
h5 = semilogy(mid(t),rads(:,5), 'Marker','x', 'color', 'k',LW,lw);

% div-con is behind bary. Plot it again on top of bary to be visible.
%h6 = semilogy(mid(t),rads(:,1),'xr', LW,lw);
h6 = semilogy(mid(t),rads(:,1), 'Marker','x', 'color', 'r',LW,lw);

xlabel('x', FS, fs)
ylabel('radius', FS, fs)
%legend(str1, 'location', 'best')
legend([h1 h2 h3 h4 h5], str1)

% div-con is behind bary. Plot it again on top of bary to be visible.
% hold on
% semilogy(mid(t),rads(:,1),LW,lw, 'Color','r')
xlabel('x', FS, fs)
ylabel('radius', FS, fs)

subplot(1,2,2)
% Keep the handle of each parts so that I can skip all but the last one
% when putting the legend!
h1 = semilogy(mid(t),rads(:,1), 'Marker','x', 'color', 'r' ,LW,lw);
% h2 = semilogy(mid(t),rads(:,1), 'color', [0 1 0], LW,lw);
hold on
h2 = semilogy(mid(t),rads(:,2), 'm<-.',LW,lw);
h3 = semilogy(mid(t),rads(:,3), 'color', [0.2 0.8 0] ,LW,lw);
h4 = semilogy(mid(t),rads(:,4), 'color', [0 1 0], LW,lw);
h5 = semilogy(mid(t),rads(:,5), 'marker', 'x', 'color', 'k', LW,lw);

% h5 = semilogy(mid(t),rads(:,5), 'oy',LW,lw);
% h6 = semilogy(mid(t),rads(:,6), 'color', [0.8 0.8 0.2], LW,lw);

h6 = semilogy(mid(t),rads(:,6), '-.', 'color', [0.85 0.33 0.1],LW,lw);
legend(h6, 'bary \cap d-cos-acos')


xlabel('x', FS, fs)
ylabel('radius', FS, fs)
text(-3.5,1.7e-6,['degree = ' num2str(size(intCoeffs,1)-1) ...
    ',  rad(c) = ' num2str(maxrad_coeffs,2) ',    l = ' num2str(l) ',   rad(x) = ' num2str(maxrad_pts,2)], FS, fs)
shg
print(gcf,'-depsc','/Users/user/Desktop/My work/git/ver-cheb-eval/draft/figures/ex9_1');

%%
figure
ax = subplot(1,2,1);
i = 1;
h=bar(i,avgDigits(i));
set(h,'FaceColor','r');
%'cos-acos'

hold on
i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor','m');
%'d-div-con'


i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', [0.2 0.8 0]);
%ICA-eig

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor','g');
%ICA-eig-err

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor','k');
%'bary'

% i = i+1;
% h=bar(i,avgDigits(i));
% set(h,'FaceColor','y');
% %'para'
% 
% i = i+1;
% h=bar(i,avgDigits(i));
% set(h,'FaceColor', [0.8 0.8 0.2]);
% %'Lohner QR'

i = i+1;
h=bar(i,avgDigits(i));
set(h,'FaceColor', [0.85 0.33 0.1]);
%'bary \cap d-div-con'

xlim([0,i+1])
ylabel('average correct digits', FS, fs)

set(ax,'xtick',1:i);
xticklabels(ax, {'d-cos-acos', 'd-div-con', 'ICA-eig', 'ICA-eig-err', 'bary', 'bary \cap d-cos-acos'});
ax.XTickLabelRotation = 45;

% Remove unnecessary parts of the second subplot
posDigInd = find( avgDigits .* (avgDigits>0) ); %index of avgDigits which are positive
m = min(avgDigits(posDigInd));
M = max(avgDigits(posDigInd));
ylim([m-0.1, M+0.1])


ax = subplot(1,2,2);
i = 1;
h=bar(i,times(i));
set(h,'FaceColor','r');
%'cos-acos'

hold on
i = i+1;
h=bar(i,times(i));
set(h,'FaceColor','m');
%'div-con'

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', [0.2 0.8 0]);
%ICA-eig

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor','g');
%ICA-eig-err

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor','k');
%'bary'

% i = i+1;
% h=bar(i,times(i));
% set(h,'FaceColor','y');
% %'para'
% 
% i = i+1;
% h=bar(i,times(i));
% set(h,'FaceColor', [0.8 0.8 0.2]);
% %'Lohner QR'

i = i+1;
h=bar(i,times(i));
set(h,'FaceColor', [0.85 0.33 0.1]);
%'bary \cap d-div-con'

xlim([0,i+1])
ylabel('time (sec)', FS, fs)
ax = gca;
set(ax,'xtick',1:i);
xticklabels(ax, {'d-cos-acos', 'd-div-con', 'ICA-eig', 'ICA-eig-err', 'bary', 'bary \cap d-cos-acos'});
ax.XTickLabelRotation = 45;
text(-9,31.5,['degree = ' num2str(size(intCoeffs,1)-1) ',  rad(c) = ' num2str(maxrad_coeffs,2) ',    l = ' num2str(l) ',   rad(x) = ' num2str(maxrad_pts,2)], FS, fs)
